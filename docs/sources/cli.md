# Command line

```sh
obdeect sources run SOURCE [options]
obdeect sources check SOURCE [options]
```

`SOURCE` is a run file (`.toml` or `.json`) or a preset name ([presets](presets.md)). A preset
takes its own options after its name; `obdeect sources run <preset> --help` lists them.

## `run`

Validates the run description, then writes the photon stream. Without `-o`, the stream goes to
standard output.

```sh
obdeect sources run run.toml -o out.pst
obdeect sources run run.toml --seed 7 --set 'events[0].repeat=5' -o out.pst
```

## `check`

Validates the run description and prints its expected values without writing photons: the
number of events and telescopes, the expected candidates and weights, and one row per emission.

```sh
obdeect sources check run.toml
obdeect sources check run.toml --json
```

`--json` prints the normalized run description instead: every default filled in. It is a valid
run file.

## Options

| Option | Command | Meaning |
| --- | --- | --- |
| `-o FILE` | `run` | photon stream output (`.pst`); default and `-`: standard output |
| `--seed N` | `run` | run seed; shortcut for `--set run.seed=N` |
| `--max-rows N` | `run` | refuse a run expected to write more than N photon rows (all telescopes, 40 bytes per row); default 1e8 |
| `--set KEY=VALUE` | both | change one value of the run description before use; repeatable |
| `--json` | `check` | print the normalized run description as JSON |

`KEY` is a path into the run description: a dotted name, and `[i]` for list items, for example
`run.seed` or `events[0].repeat`. `VALUE` is read as JSON; text that is not valid JSON is a
string:

```sh
obdeect sources run run.toml --set run.seed=7 --set run.metadata.note=test -o out.pst
```

Quote a path with brackets (`'events[0].repeat=5'`), so the shell does not expand it.

Errors go to standard error, start with `obdeect sources run:` or `obdeect sources check:`, and
exit with status 1. A message names the key that is wrong.
