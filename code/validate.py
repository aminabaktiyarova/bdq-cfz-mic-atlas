"""
Tests the pre-specified predictions on the held-out half. Runs once.

Run from the project root with the virtual environment active:

    python code/validate.py

This module reads docs/PRE_REGISTRATION.md, refuses to run if that file is not
committed, and evaluates each registered prediction against CRyPTIC-v2.0 using
the criteria fixed in that document. It changes no estimate and chooses no new
comparison.

Two guards, both deliberate.

It will not run against an uncommitted or modified pre-registration. The
document is only evidence of anything if it existed, in the state being tested
against, before the test was run. Committing it first is what makes that
checkable by someone who was not present.

It records its own outcome to outputs/validation_report.txt and leaves it there.
A failed prediction is reported as a failure. It is not explained away, and it
is not replaced with a comparison that worked.

On what this can and cannot show. The exploratory analyses that preceded the
split used the full cohort, so the held-out samples were seen in aggregate
before the comparisons were chosen. The predictions themselves were estimated on
the discovery half alone and the criteria were fixed in advance, so a successful
test shows an effect estimated in the earlier collection holding in a later and
partly different one. It does not show that the comparisons were selected blind.
docs/PRE_REGISTRATION.md states this and any write-up must repeat it.

Output: outputs/validation_report.txt
"""

import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import cohort  # noqa: E402
import discovery  # noqa: E402

PREREG = Path("docs/PRE_REGISTRATION.md")
REPORT = Path("outputs/validation_report.txt")
VALIDATION = "CRyPTIC-v2.0"
SEED = 20260101

_lines = []


def say(text=""):
    print(text)
    _lines.append(text)


def write_report():
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(_lines) + "\n")


def die(message):
    say(f"\nSTOPPED: {message}")
    write_report()
    sys.exit(1)


def git(*arguments):
    result = subprocess.run(
        ["git", *arguments], capture_output=True, text=True, check=False
    )
    return result.returncode, result.stdout.strip(), result.stderr.strip()


def check_prereg_committed():
    """Refuse to proceed unless the pre-registration is committed and unmodified."""
    if not PREREG.is_file():
        die(f"{PREREG} does not exist. Run code/discovery.py first.")

    code, _, _ = git("rev-parse", "--git-dir")
    if code != 0:
        die("this is not a git repository, so the pre-registration cannot be "
            "shown to predate the test")

    code, tracked, _ = git("ls-files", "--error-unmatch", str(PREREG))
    if code != 0:
        die(f"{PREREG} is not committed. Commit and push it before running this, "
            "because its timestamp is the only evidence it preceded the test.")

    code, diff, _ = git("diff", "HEAD", "--", str(PREREG))
    if diff:
        die(f"{PREREG} has uncommitted changes. The committed version is what the "
            "test is accountable to. Commit the changes or revert them.")

    _, commit, _ = git("log", "-1", "--format=%H %ad %s", "--date=iso", "--", str(PREREG))
    return commit


def parse_predictions(text):
    """Read the registered predictions back out of the committed document."""
    predictions = []
    current = None
    for line in text.splitlines():
        heading = re.match(r"^### (P\d+)\. (.+?), (BDQ|CFZ)\s*$", line)
        if heading:
            current = {"id": heading.group(1), "comparison": heading.group(2),
                       "drug": heading.group(3)}
            continue
        if current is None:
            continue

        binary = re.search(
            r"\*\*(P\d+a), resistance\.\*\* The held-out odds ratio against the "
            r"reference group will fall between ([\d.]+) and ([\d.]+)\.", line)
        if binary:
            predictions.append({**current, "id": binary.group(1), "kind": "binary",
                                "low": float(binary.group(2)),
                                "high": float(binary.group(3))})
            continue

        mic = re.search(
            r"\*\*(P\d+b), MIC\.\*\* The held-out fitted mean log2 MIC will "
            r"(exceed|fall below) the reference group's, by between "
            r"(-?[\d.]+) and (-?[\d.]+) doublings\.", line)
        if mic:
            predictions.append({**current, "id": mic.group(1), "kind": "mic",
                                "direction": mic.group(2),
                                "low": float(mic.group(3)),
                                "high": float(mic.group(4))})
    return predictions


def group_names(comparison):
    if comparison == "Rv0678 loss of function":
        return [f"Rv0678 {c}" for c in cohort.LOF_CLASSES]
    return [comparison]


def main():
    import numpy as np
    import pandas as pd

    rng = np.random.default_rng(SEED)

    say("=" * 72)
    say("Validation on the held-out half")
    say("=" * 72)

    commit = check_prereg_committed()
    say(f"\nPre-registration commit: {commit}")

    predictions = parse_predictions(PREREG.read_text())
    if not predictions:
        die("no predictions could be read from the pre-registration")
    say(f"Registered predictions: {len(predictions)}")
    for prediction in predictions:
        say(f"  {prediction['id']}: {prediction['comparison']}, {prediction['drug']},"
            f" {prediction['kind']}")

    df = discovery.load_cohort()
    say(f"\nHeld-out half: {int((df.DATASET == VALIDATION).sum()):,} samples")

    say("\n" + "=" * 72)
    say("Results")
    say("=" * 72)

    outcomes = []
    for prediction in predictions:
        drug = prediction["drug"]
        held_out = discovery.restrict(df, drug, VALIDATION)
        held_out = held_out[(~held_out.mmpL5_LOF) | held_out.GROUP.eq("reference")]
        groups = group_names(prediction["comparison"])

        say(f"\n{prediction['id']}. {prediction['comparison']}, {drug},"
            f" {prediction['kind']} test")
        say(f"   predicted range: {prediction['low']:.2f} to {prediction['high']:.2f}")

        if prediction["kind"] == "binary":
            result = discovery.binary_effect(held_out, groups, drug, rng)
            if not result or result["odds_ratio"] is None:
                say("   could not be estimated in the held-out half")
                outcomes.append({**prediction, "observed": None, "verdict": "not estimable"})
                continue
            observed = result["odds_ratio"]
            say(f"   held-out: {result['exposed']} samples in {result['clusters']} clusters,"
                f" {result['resistant']} resistant")
            say(f"   reference: {result['reference']:,} samples,"
                f" {result['reference_resistant']} resistant")
            say(f"   observed odds ratio: {observed:.1f}"
                f"  (its own 95% interval {result['ci_low']:.1f} to {result['ci_high']:.1f})")
            inside = prediction["low"] <= observed <= prediction["high"]
            direction_ok = observed > 1
        else:
            result = discovery.mic_shift(held_out, groups, drug, rng)
            if not result:
                say("   could not be estimated in the held-out half")
                outcomes.append({**prediction, "observed": None, "verdict": "not estimable"})
                continue
            observed = result["shift"]
            say(f"   held-out: {result['exposed']} samples in {result['clusters']} clusters")
            say(f"   reference mean {result['reference_mean']:.2f},"
                f" group mean {result['exposed_mean']:.2f}")
            say(f"   observed shift: {observed:.2f} doublings"
                f"  (its own 95% interval {result['ci_low']:.2f} to {result['ci_high']:.2f})")
            inside = prediction["low"] <= observed <= prediction["high"]
            direction_ok = (observed > 0 if prediction["direction"] == "exceed"
                            else observed < 0)

        verdict = "SUPPORTED" if (inside and direction_ok) else "FAILED"
        reason = ""
        if not direction_ok:
            reason = "wrong direction"
        elif not inside:
            reason = ("below the predicted range" if observed < prediction["low"]
                      else "above the predicted range")
        say(f"   verdict: {verdict}" + (f"  ({reason})" if reason else ""))
        outcomes.append({**prediction, "observed": float(observed), "verdict": verdict})

    say("\n" + "=" * 72)
    say("Summary")
    say("=" * 72)
    table = pd.DataFrame([{
        "prediction": o["id"], "comparison": o["comparison"], "drug": o["drug"],
        "test": o["kind"],
        "predicted": f"{o['low']:.2f} to {o['high']:.2f}",
        "observed": "" if o["observed"] is None else round(o["observed"], 2),
        "verdict": o["verdict"],
    } for o in outcomes])
    say(table.to_string(index=False))

    supported = sum(1 for o in outcomes if o["verdict"] == "SUPPORTED")
    say(f"\n{supported} of {len(outcomes)} predictions supported.")
    say("\nThis result stands as recorded. No estimate in the pre-registration is")
    say("revised in light of it, and no comparison outside the registered set is")
    say("substituted for one that failed.")

    write_report()
    say(f"\nWritten to {REPORT}")


if __name__ == "__main__":
    main()
