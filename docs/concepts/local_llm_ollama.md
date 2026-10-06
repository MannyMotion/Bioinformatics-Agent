# Local LLM (Ollama)

**What it is:** Ollama runs open-weight language models on your own computer. We use `llama3.2:3b`, a 3-billion-parameter model, called via the `ollama run` command.

**Why local:** Free (no API cost), private (genomic and patient-type data never leaves the machine), works offline, no API key to leak.

**Trade-offs (be honest about these):**
- A 3B model is much weaker at reasoning and biology than large hosted models. It can be repetitive, or confidently wrong.
- Slow: roughly 10–50 s per answer on a laptop CPU.
- Can't run on free web hosting, so a public demo needs pre-written explanations (Lesson L4).
- We call the CLI through `subprocess` and strip ANSI escape codes from its output. An HTTP client would be cleaner and allow streaming.
- A 500-character cap on user questions limits prompt abuse, but isn't real prompt-injection defence.
- Anything the LLM writes about your results is an *interpretation*. The statistics come from the pipeline, not the model.

**Where in our code:** `_call_ollama`, `explain_results`, `answer_question` in [explainer.py](../../src/bioagent/agent/explainer.py). If Ollama isn't available, we fall back to the pipeline's template text.

**Likely interview questions**
1. *Why a local LLM instead of an API?* — Cost and data privacy: bioinformatics inputs can be sensitive. I accept lower quality for that. For production I'd offer a configurable provider with proper data agreements.
2. *How do you know the LLM isn't making things up?* — I don't fully, and I say so. It only explains; the numbers come from deterministic code. I show the pipeline output next to the text and I'm working on feeding it actual result tables so claims can be checked.
