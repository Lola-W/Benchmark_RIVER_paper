#!/usr/bin/env python3
"""Create path-resolved working copies without running or submitting analyses.

Source paths beginning ./work/, ./inputs/, or ./resources/ are relative to
this package. The work tree retains the source directory layout so scripts
that change directory and submit dependent jobs keep their expected paths.
Existing files are never overwritten; a repeated identical preparation is OK.
"""
import argparse
from pathlib import Path
import re
import sys


def prepare(root):
    root = root.resolve()
    # These legacy shell snippets interpolate paths without uniform quoting.
    if re.search(r"[^A-Za-z0-9_./-]", str(root)):
        raise ValueError("Use a package location containing only letters, digits, _, -, ., and /.")
    updates = []
    for source in sorted((root / 'source').rglob('*')):
        if not source.is_file():
            continue
        target = root / 'work' / source.relative_to(root / 'source')
        content = source.read_bytes().decode()
        content = re.sub(r'\./(work|inputs|resources)(?=/|[\"\'\s}]|$)',
                         lambda m: str(root / m.group(1)), content)
        if target.exists() and target.read_bytes() != content.encode():
            raise ValueError(f"Refusing to overwrite existing working file: {target}")
        updates.append((source, target, content))
    # Check all conflicts before creating files.
    for source, target, content in updates:
        target.parent.mkdir(parents=True, exist_ok=True)
        if not target.exists():
            target.write_bytes(content.encode())
            target.chmod(source.stat().st_mode)
        # Slurm opens these log paths before the script can create directories.
        for match in re.finditer(r'^#SBATCH\s+--(?:output|error)=(\S+)', content, re.M):
            log = Path(match.group(1))
            if log.is_absolute() and root / 'work' in log.parents:
                log.parent.mkdir(parents=True, exist_ok=True)
    return len(updates)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    try:
        count = prepare(Path(__file__).resolve().parents[1])
    except (ValueError, OSError) as exc:
        sys.exit(str(exc))
    print(f"Prepared {count} working files. No analyses or jobs were run.")


if __name__ == '__main__':
    main()
