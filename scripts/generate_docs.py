"""
Generate docs by concatenating prompt blocks and (optionally) extracted code docstrings.

This offline-friendly script does not call external LLMs. It provides a simple
mechanism to keep docs in sync with prompt blocks and can be replaced with a
Codex CLI integration in environments where that is available.
"""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def _read(path: Path) -> str:
    return path.read_text() if path.exists() else ""


def main() -> None:
    prompts_dir = ROOT / "prompts" / "doc_blocks"
    docs_dir = ROOT / "docs"
    docs_dir.mkdir(parents=True, exist_ok=True)

    # Simple concatenation refresh (placeholder for Codex CLI generation)
    mapping = {
        "Documentation.md": [
            "01_mcp_role.md.gpt",
            "02_target_prediction_nondocking.md.gpt",
            "06_provenance.md.gpt",
        ],
        "Methods.md": [
            "03_comparators_geminimol.md.gpt",
            "05_pharmacophore.md.gpt",
        ],
        "Data-Sources.md": [],
        "Validation.md": [],
    }

    for doc, blocks in mapping.items():
        out = docs_dir / doc
        if blocks:
            text = "\n\n".join(_read(prompts_dir / b) for b in blocks)
            out.write_text(text)
        else:
            # Keep existing content if present
            if not out.exists():
                out.write_text("")

    print("Docs refreshed in docs/")


if __name__ == "__main__":
    main()
