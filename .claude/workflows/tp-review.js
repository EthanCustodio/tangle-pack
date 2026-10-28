export const meta = {
  name: 'tp-review',
  description: 'Review a tanglepack diff through four lenses and keep only findings two skeptics fail to refute',
  whenToUse: 'Called by /review, and by tp-implement after the change is made. args: {base, worktree}',
  phases: [
    { title: 'Find', detail: 'one reader per lens over the diff' },
    { title: 'Verify', detail: 'two refuters per finding' },
  ],
}

const { base, worktree } = args
const WHERE = `The code under review is the git worktree at ${worktree}; run every command there ` +
  `(cd ${worktree} first). The change is \`git diff ${base}\` (committed and uncommitted). ` +
  `Read whole functions around each hunk, not only the hunk.`

const LENSES = [
  { key: 'physics', prompt: 'Check the change against the physical laws and invariants in CLAUDE.md: only unstable x stable crossings, no same-stability pair treated as real, area preservation, holes propagate backward only, the topology never calls the map outside the listed exceptions, scaling only through Trellis.image_cdist, no absolute cdist slack in topology, asserted invariants kept, the numerics/topology -> loom dependency rule.' },
  { key: 'correctness', prompt: 'Look for correctness bugs: wrong logic, off-by-one, wrong branch key or orbit index, ordering assumptions, None/empty cases the invariants do not rule out, stale caches across generations, behaviour silently changed by a refactor.' },
  { key: 'style', prompt: 'Check the change against .claude/skills/house-style/SKILL.md and the Coding Style section of CLAUDE.md: needless helpers or wrappers, defensive code for impossible states, narration comments, names from plumbing rather than the paper, missing docstrings on public API, leftover debug output. Report only what a careful author would actually change.' },
  { key: 'tests', prompt: 'Check the tests: is new or changed behaviour exercised, do new tests have teeth (would they fail if the change were reverted?), do they follow the test policy in tests/conftest.py (no pinned case-specific words, letters, matrices or counts; laws over cases.LAW_CASES), and did any test get weakened or deleted without reason.' },
]

const FINDINGS = {
  type: 'object',
  properties: {
    findings: {
      type: 'array',
      items: {
        type: 'object',
        properties: {
          file: { type: 'string' },
          line: { type: 'integer' },
          severity: { type: 'string', enum: ['high', 'medium', 'low'] },
          summary: { type: 'string' },
          failure_scenario: { type: 'string', description: 'concrete inputs/state -> wrong result, or the concrete reader cost for style' },
        },
        required: ['file', 'line', 'severity', 'summary', 'failure_scenario'],
      },
    },
  },
  required: ['findings'],
}

const VERDICT = {
  type: 'object',
  properties: { refuted: { type: 'boolean' }, reason: { type: 'string' } },
  required: ['refuted', 'reason'],
}

const RANK = { high: 0, medium: 1, low: 2 }

const results = await pipeline(
  LENSES,
  lens => agent(`${WHERE}\n\nLens: ${lens.key}. ${lens.prompt}\n\nReturn only real findings, each anchored to a file and line of the NEW code. An empty list is a fine answer.`,
    { label: `find:${lens.key}`, phase: 'Find', schema: FINDINGS }),
  (found, lens) => parallel((found?.findings ?? []).map(f => () =>
    parallel([0, 1].map(i => () => agent(
      `${WHERE}\n\nA reviewer (${lens.key} lens) claims:\n${f.file}:${f.line} [${f.severity}] ${f.summary}\nScenario: ${f.failure_scenario}\n\n` +
      `Try to REFUTE it (skeptic ${i + 1}): read the code, check whether the scenario can occur given CLAUDE.md's invariants, run a quick probe if that settles it. ` +
      `refuted=true if the claim is wrong, already handled, or out of the diff's scope.`,
      { label: `verify:${f.file.split('/').pop()}:${f.line}`, phase: 'Verify', schema: VERDICT, effort: 'medium' })))
      .then(votes => ({ ...f, lens: lens.key, survives: votes.filter(Boolean).length === 2 && votes.every(v => v && !v.refuted) }))
  )),
)

const all = results.filter(Boolean).flat().filter(Boolean)
const confirmed = all.filter(f => f.survives).sort((a, b) => RANK[a.severity] - RANK[b.severity])
log(`${all.length} findings, ${confirmed.length} survived verification`)
return { confirmed, refuted_count: all.length - confirmed.length }
