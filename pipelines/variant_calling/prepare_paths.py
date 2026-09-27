#!/usr/bin/env python3
"""Resolve relative YAML/TSV paths into disposable .runtime/ input files.

Only values beginning './' or '../' are resolved. Non-path values, caller
parameters, PATH commands, and container-internal paths are unchanged.
"""
import argparse
from pathlib import Path


def resolve(value, directory):
    if isinstance(value, dict):
        return {key: resolve(item, directory) for key, item in value.items()}
    if isinstance(value, list):
        return [resolve(item, directory) for item in value]
    if isinstance(value, str) and value.startswith(("./", "../")):
        # Use an absolute lexical path, preserving user-created symlinks.
        import os
        return os.path.abspath(directory / value)
    return value


def prepare(source):
    source = Path(source).absolute()
    target = source.parent / ".runtime" / source.name
    if source.suffix == ".tsv":
        lines = []
        for line in source.read_text().splitlines(keepends=True):
            if line.startswith("#") or not line.strip():
                lines.append(line)
            else:
                fields = line.rstrip("\r\n").split("\t")
                lines.append("\t".join(resolve(field, source.parent) for field in fields) + "\n")
        content = "".join(lines)
    elif source.suffix in (".yaml", ".yml"):
        import yaml
        config = resolve(yaml.safe_load(source.read_text()), source.parent)
        # These sample sheets also contain paths relative to their own directory.
        if isinstance(config, dict):
            for key in ("input_files", "samples_tsv", "input_paths"):
                if key in config:
                    config[key] = str(prepare(config[key]))
        content = yaml.safe_dump(config, sort_keys=False)
    else:
        raise ValueError("Expected a YAML configuration or TSV sample sheet")
    target.parent.mkdir(exist_ok=True)
    target.write_text(content)
    return target


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", help="Relative YAML configuration or TSV sample sheet")
    args = parser.parse_args()
    print(prepare(args.source))
