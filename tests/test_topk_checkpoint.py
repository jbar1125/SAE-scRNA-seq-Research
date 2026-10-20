"""CPU checks for opt-in TopK training checkpoints; no research data required.

Run: python3 tests/test_topk_checkpoint.py
"""
import contextlib
import hashlib
import io
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import numpy as np
import torch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from causal_pipeline import _topk_input_metadata, train_topk_sae


def _data():
    return np.random.default_rng(19).normal(size=(73, 11)).astype(np.float32)


def _train(rep, **kwargs):
    params = {"latent_dim": 17, "k": 4, "batch": 16, "lr": 0.003,
              "seed": 9, "device": "cpu"}
    params.update(kwargs)
    return train_topk_sae(rep, **params)


def _equal(left, right):
    if isinstance(left, torch.Tensor):
        assert torch.equal(left, right), "Tensor state differs"
    elif isinstance(left, dict):
        assert left.keys() == right.keys()
        for key in left:
            _equal(left[key], right[key])
    elif isinstance(left, (list, tuple)):
        assert len(left) == len(right)
        for a, b in zip(left, right):
            _equal(a, b)
    else:
        assert left == right, (left, right)


def _raises(error, call):
    try:
        call()
    except error:
        return
    raise AssertionError(f"Expected {error.__name__}")


def test_exact_resume_and_inference():
    rep = _data()
    with tempfile.TemporaryDirectory() as directory:
        continuous = Path(directory) / "continuous.pt"
        resumed = Path(directory) / "resumed.pt"
        encode_full, full_info = _train(rep, epochs=5, checkpoint_path=continuous)
        _train(rep, epochs=2, checkpoint_path=resumed)
        encode_resumed, resumed_info = _train(rep, epochs=5, checkpoint_path=resumed, resume=True)
        full = torch.load(continuous, weights_only=True)
        split = torch.load(resumed, weights_only=True)
        for key in ("enc", "dec", "optimizer", "mu", "rng_cpu", "rng_cuda", "rng_numpy",
                    "input", "training_params", "completed_epochs", "target_epochs"):
            _equal(full[key], split[key])
        np.testing.assert_array_equal(encode_full(rep), encode_resumed(rep))
        assert full_info["epochs_completed"] == resumed_info["epochs_completed"] == 5
        assert resumed_info["checkpoint_path"] == str(resumed.resolve())

        # A completed checkpoint must reconstruct inference without any optimizer steps.
        encode_loaded, info = _train(rep, epochs=5, checkpoint_path=resumed, resume=True)
        np.testing.assert_array_equal(encode_full(rep), encode_loaded(rep))
        assert info["epochs_completed"] == 5

        # Enabling checkpoints must not change the default training trajectory or output.
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            encode_default, default_info = _train(rep, epochs=5)
        assert captured.getvalue() == ""
        assert set(default_info) == {"latent_dim", "k", "d", "device"}
        np.testing.assert_array_equal(encode_default(rep), encode_full(rep))
    print("PASS exact uninterrupted vs 2+3 resumed weights, optimizer, RNG and inference")


def test_refuse_overwrite_and_mismatches():
    rep = _data()
    with tempfile.TemporaryDirectory() as directory:
        checkpoint = Path(directory) / "sae.pt"
        _train(rep, epochs=2, checkpoint_path=checkpoint)
        original_bytes = checkpoint.read_bytes()
        _raises(FileExistsError, lambda: _train(rep, epochs=3, checkpoint_path=checkpoint))
        _raises(FileNotFoundError, lambda: _train(rep, epochs=3, resume=True,
                                                 checkpoint_path=Path(directory) / "missing.pt"))
        _raises(ValueError, lambda: _train(rep, epochs=3, resume=True))
        _raises(ValueError, lambda: _train(rep, epochs=1, checkpoint_path=checkpoint, resume=True))
        for override in ({"latent_dim": 18}, {"k": 3}, {"batch": 17}, {"lr": 0.001}, {"seed": 8}):
            _raises(ValueError, lambda override=override: _train(
                rep, epochs=3, checkpoint_path=checkpoint, resume=True, **override))
        changed = rep.copy()
        changed[0, 0] += np.float32(0.001)
        for changed_input in (changed, rep[:-1], rep.astype(np.float64), rep[::-1]):
            _raises(ValueError, lambda changed_input=changed_input: _train(
                changed_input, epochs=3, checkpoint_path=checkpoint, resume=True))
        assert checkpoint.read_bytes() == original_bytes
        _raises(ValueError, lambda: _train(rep, epochs=1, checkpoint_every=0,
                                           checkpoint_path=Path(directory) / "invalid.pt"))
    print("PASS explicit-resume requirement, existing-file protection, input/parameter refusals")


def test_atomic_save_and_checkpoint_schedule():
    rep = _data()
    with tempfile.TemporaryDirectory() as directory:
        checkpoint = Path(directory) / "sae.pt"
        captured = io.StringIO()
        with contextlib.redirect_stdout(captured):
            _train(rep, epochs=11, checkpoint_path=checkpoint)
        lines = captured.getvalue().splitlines()
        assert len(lines) == 2 and "epoch 10/11" in lines[0] and "epoch 11/11" in lines[1]
        assert torch.load(checkpoint, weights_only=True)["completed_epochs"] == 11
        original_bytes = checkpoint.read_bytes()
        with patch("torch.save", side_effect=OSError("simulated write failure")):
            _raises(OSError, lambda: _train(rep, epochs=11, checkpoint_path=checkpoint, resume=True))
        assert checkpoint.read_bytes() == original_bytes, "Failed atomic save damaged previous checkpoint"
        assert list(Path(directory).glob("*.tmp")) == [], "Failed save left a temporary file"
    print("PASS ten-epoch/final checkpoint schedule and atomic-save failure protection")


def test_streaming_input_hash():
    # Noncontiguous input larger than one chunk verifies logical array ordering and
    # bounded memoryview updates, without imposing a large extra bytes allocation.
    rep = np.arange(5000 * 600, dtype=np.float32).reshape(5000, 600)[:, ::2]
    expected = hashlib.sha256(rep.tobytes(order="C")).hexdigest()
    original_sha256 = hashlib.sha256
    update_sizes = []

    class RecordedDigest:
        def __init__(self):
            self.digest = original_sha256()

        def update(self, chunk):
            assert isinstance(chunk, memoryview)
            update_sizes.append(chunk.nbytes)
            self.digest.update(chunk)

        def hexdigest(self):
            return self.digest.hexdigest()

    with patch("hashlib.sha256", RecordedDigest):
        actual = _topk_input_metadata(rep)
    assert actual == {"shape": [5000, 300], "dtype": rep.dtype.str, "sha256": expected}
    assert len(update_sizes) > 1 and max(update_sizes) <= 4 * 1024 * 1024
    assert sum(update_sizes) == rep.nbytes
    print("PASS chunked SHA256 for noncontiguous float32 input")


if __name__ == "__main__":
    torch.set_num_threads(1)
    test_exact_resume_and_inference()
    test_refuse_overwrite_and_mismatches()
    test_atomic_save_and_checkpoint_schedule()
    test_streaming_input_hash()
    print("ALL TOPK CHECKPOINT TESTS PASSED")
