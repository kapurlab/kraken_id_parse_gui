"""Saved target names live outside the checkout, and keep people's notes.

The "Target taxon" dropdown is the shipped config/taxa.yaml plus the names
users save. Saving used to append to config/taxa.yaml, a tracked file: the
checkout went dirty, `bdtools update` / `sync` refused to move the tool, and on
a site install (checkout owned by the admin) the save failed. Saved names now go
to <config dir>/kraken_id_parse_gui/target_names.yaml, and names an older
release appended are carried there once.

Run directly:  <env>/bin/python backend/app/test_target_names.py
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import target_names as T  # noqa: E402

FAILED = 0


def check(label, cond):
    global FAILED
    if cond:
        print(f"  OK  {label}")
    else:
        FAILED += 1
        print(f"  FAIL {label}")


SHIPPED = (
    "# Kraken ID Parse — taxon search names\n"
    "\n"
    "- Mycobacterium tuberculosis complex\n"
    "- Orbivirus     # BTV/EHD, beta\n"
)


def _repo(tmp: Path) -> Path:
    root = tmp / "kraken_id_parse_gui"
    (root / "config").mkdir(parents=True)
    (root / "config" / "taxa.yaml").write_text(SHIPPED, encoding="utf-8")
    git = ["git", "-C", str(root)]
    subprocess.run(git + ["init", "-q"], check=True)
    subprocess.run(git + ["add", "config/taxa.yaml"], check=True)
    subprocess.run(git + ["-c", "user.name=t", "-c", "user.email=t@example.invalid",
                          "commit", "-q", "-m", "ship"], check=True)
    return root


def _status(root: Path) -> str:
    return subprocess.run(["git", "-C", str(root), "status", "--porcelain"],
                          capture_output=True, text=True).stdout.strip()


def main() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="target_names_"))
    os.environ["XDG_CONFIG_HOME"] = str(tmp / "xdg")
    root = _repo(tmp)
    saved = T.saved_path()
    check("saved list lives in the per-user config folder",
          saved == tmp / "xdg" / "kraken_id_parse_gui" / "target_names.yaml")

    print("list")
    check("shipped names come first, as written",
          T.names(root) == ["Mycobacterium tuberculosis complex", "Orbivirus"])

    print("add")
    r = T.add(root, "  Brucella   abortus ")
    check("a new name is saved, whitespace tidied", r == {"name": "Brucella abortus", "added": True})
    check("...to the saved file, not the checkout", saved.is_file() and _status(root) == "")
    check("...after the shipped names", T.names(root)[-1] == "Brucella abortus")
    check("entries say where each name comes from",
          [e["source"] for e in T.entries(root)] == ["shipped", "shipped", "saved"])
    check("a duplicate (any case) is not added twice",
          T.add(root, "orbivirus") == {"name": "Orbivirus", "added": False}
          and T.add(root, "BRUCELLA ABORTUS")["added"] is False)
    try:
        T.add(root, "   ")
        check("an empty name is refused", False)
    except ValueError:
        check("an empty name is refused", True)

    # A note kept by hand in the saved file survives later adds and removes.
    saved.write_text(saved.read_text(encoding="utf-8")
                     + "# for the 2026 outbreak — remove after the season\n", encoding="utf-8")
    T.add(root, "Weird: name")
    T.add(root, "Salmonella enterica")
    text = saved.read_text(encoding="utf-8")
    check("a hand-written note survives an add",
          "# for the 2026 outbreak — remove after the season" in text)
    check("a name needing quotes is quoted", '- "Weird: name"' in text)
    check("...and reads back as itself", "Weird: name" in T.names(root))

    print("remove")
    check("a saved name can be removed", T.remove(root, "weird: NAME") == {"removed": True})
    text = saved.read_text(encoding="utf-8")
    check("...only that line goes", "Weird" not in text and "Salmonella enterica" in text
          and "# for the 2026 outbreak" in text and text.startswith("# Kraken ID Parse"))
    try:
        T.remove(root, "Orbivirus")
        check("a shipped name cannot be removed", False)
    except ValueError:
        check("a shipped name cannot be removed", True)
    check("removing a name that is not there is a no-op",
          T.remove(root, "Nothing here") == {"removed": False})

    print("carry-over from older releases")
    shipped = root / "config" / "taxa.yaml"
    shipped.write_text(SHIPPED + "- Poxviridae\n- Salmonella enterica\n", encoding="utf-8")
    check("(an older release left the checkout dirty)", _status(root) != "")
    res = T.migrate_checkout_additions(root)
    check("appended names are carried into the saved list",
          res["carried"] == ["Poxviridae"] and "Poxviridae" in T.names(root))
    check("...a name already saved is not added twice",
          T.names(root).count("Salmonella enterica") == 1)
    check("the shipped file is back to the committed text", shipped.read_text(encoding="utf-8") == SHIPPED)
    check("...so the checkout is clean and updates can move it", _status(root) == "" and res["restored"])
    check("running it again does nothing", T.migrate_checkout_additions(root) == {"carried": [], "restored": False})

    edited = SHIPPED.replace("- Orbivirus", "- Orbivirus\n- Rhabdoviridae")
    shipped.write_text(edited, encoding="utf-8")
    res = T.migrate_checkout_additions(root)
    check("an edit that is not a plain append is carried but left in place",
          res["carried"] == ["Rhabdoviridae"] and not res["restored"]
          and shipped.read_text(encoding="utf-8") == edited)

    print("Kraken report names")
    rep = tmp / "s_reportkraken.txt"
    rep.write_text(
        "  2.38\t100\t100\tU\t0\tunclassified\n"
        " 97.62\t4000\t12\tR\t1\troot\n"
        " 96.00\t3900\t20\tG\t1763\t              Mycobacterium\n"
        " 95.00\t3800\t3800\tG1\t77643\t                Mycobacterium tuberculosis complex\n"
        "  0.01\t1\t1\tS\t9606\t      Homo sapiens\n",
        encoding="utf-8")
    rows = T.report_names(rep, 0.05)
    check("names come back un-indented, root/unclassified skipped, below-floor dropped",
          [r["name"] for r in rows] == ["Mycobacterium", "Mycobacterium tuberculosis complex"])
    check("...with %, reads and rank", rows[1] == {"name": "Mycobacterium tuberculosis complex",
                                                   "percent": 95.0, "reads": 3800, "rank": "G1"})

    print("PASS" if not FAILED else f"{FAILED} FAILURE(S)")
    return 1 if FAILED else 0


if __name__ == "__main__":
    raise SystemExit(main())
