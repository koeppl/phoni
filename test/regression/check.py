#!/usr/bin/env python3
"""Regression test for the matching-statistics query initialization.

Builds an index for reference.fa in a temporary directory, queries
queries.fa, and checks every MS length against a brute-force substring
search and every pointer against the reference text.

Usage: python3 test/regression/check.py <build-dir>
"""
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent


def read_fasta(path):
    records = {}
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            name = line[1:].split()[0]
            records[name] = ''
        elif line.strip():
            records[name] += line.strip()
    return records


def read_output(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith('>'):
            rows.append([])
        elif line.strip():
            rows[-1].extend(map(int, line.split()))
    return rows


def matching_statistics(query, text):
    lengths = []
    for i in range(len(query)):
        length = 0
        while i + length < len(query) and query[i:i + length + 1] in text:
            length += 1
        lengths.append(length)
    return lengths


def main():
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    build = Path(sys.argv[1]).resolve()
    text = ''.join(read_fasta(HERE / 'reference.fa').values())
    queries = read_fasta(HERE / 'queries.fa')
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        reference = tmp / 'reference.txt'
        reference.write_text(text)
        patterns = tmp / 'queries.fa'
        patterns.write_text((HERE / 'queries.fa').read_text())
        subprocess.run([sys.executable, str(build / 'phoni'), 'build', '-r', str(reference),
                        '-w', '4', '-p', '10'], check=True, stdout=subprocess.DEVNULL)
        subprocess.run([str(build / 'test' / 'src' / 'phoni'), str(reference), '-p', str(patterns),
                        '-g', 'naive'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        lengths = read_output(Path(str(patterns) + '.lengths'))
        pointers = read_output(Path(str(patterns) + '.pointers'))
    failures = 0
    for (name, query), observed, refs in zip(queries.items(), lengths, pointers):
        expected = matching_statistics(query, text)
        bad_pointer = any(text[p:p + n] != query[i:i + n] for i, (p, n) in enumerate(zip(refs, observed)))
        if observed != expected or bad_pointer:
            failures += 1
            print(f'FAIL {name}: lengths {observed}, expected {expected}, pointers {refs}')
    print(f'{len(queries) - failures}/{len(queries)} queries passed')
    return 1 if failures else 0


if __name__ == '__main__':
    sys.exit(main())
