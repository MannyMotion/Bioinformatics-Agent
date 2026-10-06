# START HERE

Welcome back, Manny. This is the one page to open first.

## The idea in one paragraph

BioAgent already works: you upload a FASTA, a count matrix or a VCF, it works out what the file is, runs the right analysis and explains the result. The **new job** is not adding features. It is making you the person who can explain every part of it, prove it is correct, and show it off. In an interview, "I built this and here's how I know it's right" beats "look how many features it has."

## What a lesson looks like

1. **Goal.** One clear thing, stated up front.
2. **Concept + biology** in plain English.
3. **You choose the help level:** hint, step-by-step, or exact lines with each one explained.
4. **You type it and run it.** I don't paste core code. You learn by typing.
5. **I review it like a senior colleague:** what's good, what to fix, and why.
6. **You explain it back.** If it's shaky, I teach it a different way.
7. **Commit**, then the LEARNING_LOG row, the DECISIONS entry and the concept note get updated.

Say **"Lesson W1"** to begin.

## The map

| File | What it's for |
|---|---|
| [CLAUDE.md](CLAUDE.md) | The rules I follow every session (teaching mode, the validation rule) |
| [PROJECT_MEMORY.md](PROJECT_MEMORY.md) | Where the project is right now |
| [ROADMAP.md](ROADMAP.md) | The lessons, in order. Phase 0 understand, Phase 1 prove, Phase 2 show, Phase 3 grow |
| [LEARNING_LOG.md](LEARNING_LOG.md) | Your skills checklist. Honest yes / shaky / no |
| [DECISIONS.md](DECISIONS.md) | Why things are built the way they are, with interview-ready wording |
| [IDEAS.md](IDEAS.md) | Backlog of what to build later, ranked by hiring value |
| [DATA_SOURCES.md](DATA_SOURCES.md) | Real public data to validate against |
| [docs/concepts/](docs/concepts/) | One short note per method: CPM, log2FC, t-test, FDR, volcano, GC, VCF, ClinVar, RAG, Ollama |
| [docs/interview/INTERVIEW_BANK.md](docs/interview/INTERVIEW_BANK.md) | Likely questions with answers in your voice. Unearned ones are marked |
| [.claude/skills/tutor/SKILL.md](.claude/skills/tutor/SKILL.md) | How I teach |
| [.claude/skills/bioagent-map/SKILL.md](.claude/skills/bioagent-map/SKILL.md) | Where everything is in the code |

## Honest state of play

- It runs and all 23 tests pass, but those tests mostly check "does it run", not "is it right".
- The weak spots are listed in [PROJECT_MEMORY.md](PROJECT_MEMORY.md) section 5. The big one: RNA-seq is a t-test on **raw counts** with no multiple-testing correction.
- The fix is **Phase 1**. After that you can say "validated against pydeseq2", and that is the sentence that lands interviews.

## Run it

See the Environment section of [CLAUDE.md](CLAUDE.md). The short version: start uvicorn on port 8000, open `http://localhost:8000/frontend/index.html`.
