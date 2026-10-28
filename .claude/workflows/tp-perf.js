export const meta = {
  name: 'tp-perf',
  description: 'Try tanglepack optimization candidates in parallel worktrees, then benchmark base and each candidate one at a time',
  whenToUse: 'Called by /perf after plan approval. args: {target, benchmark, candidates: [{name, idea, files}], base, repo}',
  phases: [
    { title: 'Implement', detail: 'one worktree per candidate' },
    { title: 'Measure', detail: 'serial, so timings do not compete for CPU' },
  ],
}

const { target, benchmark, candidates, base, repo } = args
const PY = `MPLBACKEND=Agg ${repo}/env/bin/python`

const BUILT = {
  type: 'object',
  properties: {
    path: { type: 'string', description: 'absolute path of the worktree you worked in (pwd)' },
    branch: { type: 'string' },
    suite: { type: 'string', description: 'pytest summary line' },
    equal_results: { type: 'string', description: 'how you checked the results are unchanged, and the verdict' },
    notes: { type: 'string' },
  },
  required: ['path', 'branch', 'suite', 'equal_results', 'notes'],
}

const MEASURED = {
  type: 'object',
  properties: {
    seconds: { type: 'number' },
    peak_mib: { type: 'number' },
    raw: { type: 'string', description: 'the benchmark output lines the numbers came from' },
  },
  required: ['seconds', 'peak_mib', 'raw'],
}

phase('Implement')
const built = await parallel(candidates.map(c => () => agent(
  `Optimize tanglepack for: ${target}. Your candidate: "${c.name}": ${c.idea} (files: ${c.files.join(', ')}).\n` +
  `You are in a fresh git worktree. First \`git checkout -b perf/${c.name} ${base}\`. Python is \`${PY}\` (pytest picks up this worktree's src; scripts need PYTHONPATH=src). ` +
  `Follow .claude/skills/house-style/SKILL.md. The optimization must not change what is computed: run the full suite, and compare the benchmark's outputs (registry, crossings, words) against base where the benchmark prints them. ` +
  `Do NOT time anything; timing happens later, serially. Commit your change (end the message with "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>").`,
  { label: `build:${c.name}`, phase: 'Implement', isolation: 'worktree', schema: BUILT }
).then(b => b && { ...b, name: c.name })))

const ok = built.filter(Boolean)
log(`${ok.length}/${candidates.length} candidates built`)

// Serial on purpose: concurrent benchmarks would distort each other's timings.
phase('Measure')
const sites = [{ name: 'base', path: repo, ref: base }, ...ok.map(b => ({ name: b.name, path: b.path }))]
const table = []
for (const s of sites) {
  const where = s.ref
    ? `Create a detached worktree: \`git -C ${repo} worktree add --detach ${repo}-perf-base ${s.ref}\`, cd there, and remove it afterwards with \`git -C ${repo} worktree remove ${repo}-perf-base\`.`
    : `cd ${s.path}.`
  const m = await agent(
    `${where} Run the benchmark three times and report the median wall time and peak memory: \`${benchmark}\` (Python is \`${PY}\`; prefix PYTHONPATH=src for a script). Change nothing.`,
    { label: `measure:${s.name}`, phase: 'Measure', effort: 'low', schema: MEASURED })
  table.push({ name: s.name, ...(m ?? { seconds: null, peak_mib: null, raw: 'measurement failed' }),
    ...(ok.find(b => b.name === s.name) ?? {}) })
}

return { table }
