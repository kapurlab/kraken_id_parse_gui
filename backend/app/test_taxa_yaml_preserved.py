"""Adding a taxon must not rewrite config/taxa.yaml over its comments.

`taxa.yaml` says entries may be added "by hand", and both this GUI and the
vSNP GUI write the same file. Adding one from either used to rebuild the list
from the names the parser understood, which deleted every comment in it — the
note saying why a taxon is there and when it can come out. An add is now an
append.

Run directly:  <conda>/bin/python backend/app/test_taxa_yaml_preserved.py
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import main as M  # noqa: E402

FAILED = 0


def check(label, cond):
    global FAILED
    if cond:
        print(f"  OK  {label}")
    else:
        FAILED += 1
        print(f"  FAIL {label}")


def main() -> int:
    taxa = Path(tempfile.mkdtemp(prefix="taxa_")) / "config" / "taxa.yaml"
    taxa.parent.mkdir(parents=True)
    taxa.write_text(
        "# Kraken ID Parse — taxon search names\n"
        "#\n"
        "# add by hand or via the \"Add search name\" control in either GUI.\n"
        "\n"
        "- Mycobacterium tuberculosis     # MTBC0 reference set\n"
        "# added for the 2026 outbreak — remove after the season\n"
        "- Brucella abortus\n",
        encoding="utf-8",
    )
    M._TAXA_YAML = taxa

    M._append_taxon("Salmonella enterica")
    text = taxa.read_text(encoding="utf-8")
    check("the standalone comment survived",
          "# added for the 2026 outbreak — remove after the season" in text)
    check("the inline comment survived", "# MTBC0 reference set" in text)
    check("the header survived", text.startswith("# Kraken ID Parse"))
    check("the new name went on the end",
          text.rstrip().endswith("- Salmonella enterica"))
    check("every name reads back, in order",
          M._read_taxa() == ["Mycobacterium tuberculosis", "Brucella abortus",
                             "Salmonella enterica"])

    M._append_taxon("Weird: name")
    check("a name needing quotes is still quoted",
          '- "Weird: name"' in taxa.read_text(encoding="utf-8"))
    check("...and reads back as itself", M._read_taxa()[-1] == "Weird: name")

    missing = taxa.parent / "gone.yaml"
    M._TAXA_YAML = missing
    M._append_taxon("First entry")
    check("a file that does not exist yet is created with the header",
          missing.read_text(encoding="utf-8").startswith("# Kraken ID Parse")
          and M._read_taxa() == ["First entry"])

    print("PASS" if not FAILED else f"{FAILED} FAILURE(S)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
