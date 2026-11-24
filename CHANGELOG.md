# Changelog

## 1.1.2 (unreleased)

Error handling and robustness. Results for valid input are unchanged; runs that previously finished silently with missing or incomplete results now stop with an error.

### Upgrade notes
- A run succeeded only if the exit code is 0 and `log/<study name>.summary.tsv` says `status ok`. Failed runs exit with 1, 2 or 3 instead of 0; wrappers that relied on exit code 0 should check both.
- Runs can now stop where they previously continued: non-numeric configuration values, missing tools or reference files, a `refchr` that is not in the reference indexes, a failed command, write or read, no reads or alignments after a step, missing copy number coverage, or inconsistent intermediate files.
- `log/<study name>.log` is overwritten by each run instead of appended to; back it up before re-running with the same study name.
- The scripts are now `MitoSAlt.pl` and `MitoSAlt_SE.pl`; `MitoSAlt1.1.1.pl` and `MitoSAlt_SE1.1.1.pl` remain as links to them.
- The version is read from the `VERSION` file.

### Added
- Exit codes: 0 success, 1 invalid input or configuration, 2 a pipeline step failed, 3 no MT reads or alignments left after a step (an input or reference problem, not a negative result).
- Run summary `log/<study name>.summary.tsv`, one `key<TAB>value` per line. It is written as `status running` at the start and replaced on success or failure, so a summary left at `running` means the process was killed. Keys: `status` (running/ok/failed), `exit_code`, `message` (failed runs), `version`, `sample` (the study name), `mt_reads` (reads on the MT contig in the nuclear+MT alignment), `mt_alignments` (alignments in the MT BAM), `mt_mean_depth` (mean coverage of the MT contig, to judge whether a run without events had enough coverage), `split_reads` (lines in `.breakpoint`), `clusters` (lines in `.cluster`), `events_reported` (rows in `.tsv`), `warnings` (Perl warnings during the run). Keys for steps that did not run are omitted.
- The final `.tsv` is always written; a header-only file means no events passed the filters.
- Checks before running: arguments, configuration keys, STEPS values, numeric parameter values, tools, indexes and reference files for the enabled steps, and that `refchr` is a contig in the reference indexes. Missing configuration or FASTQ files are reported in the log and summary.
- Optional configuration key `javamem` (GB, default 100) for the BBMap `reformat` Java heap in enriched mode, previously fixed at 100 GB; same key as the nf-core/raredisease patch.
- `--version` option.
- Optional delplot.R arguments for the plot and TSV output directories.

### Changed
- One log per run, `log/<study name>.log`, also shown on the terminal. It starts with the version, command, host and working directory and holds the messages of the run and the output of every tool, including HISAT2 and BBMap, which previously wrote only to the file. Perl warnings are marked `WARNING:`; R warnings are printed as they occur.
- Every external command runs with `pipefail` and the pipeline stops on the first failure, instead of finishing with exit code 0. Temporary files are kept when a run fails.
- Perl errors and SIGINT/SIGTERM/SIGHUP end the run with exit code 2 and a summary.
- delplot.R is run with Rscript from the script directory; its messages go to the run log instead of `delplot.Rout`.
- delplot.R stops with a message when an R package is missing, instead of installing packages during a run.
- Output directories are created automatically.

### Fixed
- Writing and reading of `.breakpoint`, `.bed`, `.cluster`, the temporary files and the summary are checked, so a full disk or read error stops the run instead of leaving truncated files. The summary is written to a temporary file and renamed. Commands whose output is read directly (read naming detection, the MT alignment count) are checked with `pipefail` too.
- BBMap `reformat` (enriched mode) runs with `eoom`, so a Java out-of-memory error during the run exits non-zero instead of possibly leaving a truncated FASTQ.
- Read naming (paired `/1`, `/2` or unpaired) is detected from the last line of `tab/<study name>.tab.gz`; the run stops when that line is not an alignment row, instead of possibly parsing read names the wrong way.
- The copy number step (cn_mt) stops with a message when no random regions are generated from `gsize` (exit code 2), when there are no nuclear reads in them, e.g. for MT-only input (exit code 3), or when `sambamba depth window` writes no coverage rows (exit code 2). Previously `bedtools groupby` reported the missing input but exited 0, so a line was silently missing from `.cnmt.txt`.
- delplot.R stops when clusters have no reads in the `.breakpoint` file, instead of writing a header-only `.tsv` that looks like a negative result.
- A `plot/<study name>.pdf` from an earlier run is removed before plotting, so it cannot be mistaken for this run's result.
- `rmtmp = yes` removes this run's temporary files by exact name, instead of with the pattern `tmp_<study name>*`, which also removed the files of other runs in the same directory (e.g. `S1` removed those of `S10`).
- Host shell startup files given in `BASH_ENV`/`ENV` (e.g. an lmod setup passed through by Singularity) are no longer run for every command.
- `insall.packages` typo in delplot.R.
- README: script names, the `hsindex` key, the `lastdb` path, the alignment table file name and the `.tsv` column names now match the code.
