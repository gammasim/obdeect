"""Canonical hashes for compiled optical-model documents."""

import hashlib
import json
from typing import Any


def model_hash(model: dict[str, Any]) -> str:
    """Return the SHA-256 of canonical JSON without the document's ``sha256`` field."""
    content = {key: value for key, value in model.items() if key != "sha256"}
    text = json.dumps(content, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(text.encode()).hexdigest()
