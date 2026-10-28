export const meta = {
  name: 'tp-implement',
  description: 'Implement an approved tanglepack plan in its worktree, review it, fix confirmed findings, and check suite and coverage',
  whenToUse: 'Called by /refactor and /feature after plan approval. args: {mode, slug, repo, worktree, plan, base, base_cov, base_suite}',
  phases: [
    { title: 'Implement' },
    { title: 'Review' },
    { title: 'Fix' },
    { title: 'Check' },
  ],
}

const { mode, slug, repo, worktree, plan, base, base_cov, base_suite } = args
const PY = `MPLBACKEND=Agg ${repo}/env/bin/python`
const WHERE = `Work ONLY in the git worktree ${worktree} (branch ${mode}/${slug}); cd there before every command and edit files under that path only. ` +
  `Python is \`${PY}\` (pytest picks up the worktree's src; scripts need PYTHONPATH=src). ` +
  `Follow .claude/skills/house-style/SKILL.md. Never touch the main repo checkout.`

const MODE_RULES = {
  refactor: 'This is a REFACTOR: nothing computed may change. Write a characterization test only for touched code that no test covers.',
  feature: 'This is a FEATURE: write the tests the plan names for the new behaviour, following the test policy in tests/conftest.py.',
}

const REPORT = {
  type: 'object',
  properties: {
    commits: { type: 'array', items: { type: 'string' }, description: 'sha + subject, oldest first' },
    suite: { type: 'string', description: 'pytest summary line, e.g. "1130 passed, 1 xfailed"' },
    notes: { type: 'string', description: 'anything the author must know: deviations from the plan, scope questions, output changes' },
  },
  required: ['commits', 'suite', 'notes'],
}

phase('Implement')
const impl = await agent(
  `${WHERE}\n\n${MODE_RULES[mode]}\n\nRead the approved plan at ${plan} and carry it out exactly. The plan's file list is the scope: if the change needs another file, do the minimum and say so in notes. ` +
  `Commit in small logical commits (message style: "<Area>: <what changed>", e.g. git log --oneline in the repo; end each message with "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"). ` +
  `Run the full suite at the end (baseline was: ${base_suite}); it must not regress.`,
  { label: `implement:${slug}`, phase: 'Implement', schema: REPORT })
if (!impl) return { error: 'implementer did not finish' }

const leftovers = []
let findings = []
for (let round = 1; round <= 2; round++) {
  phase('Review')
  const review = await workflow('tp-review', { base, worktree })
  findings = review?.confirmed ?? []
  if (!findings.length) break
  phase('Fix')
  const fix = await agent(
    `${WHERE}\n\n${MODE_RULES[mode]}\n\nApply these verified review findings (round ${round}). Fix each one that is in the plan's scope (${plan}); for any you decline, say why in notes. Run the suite, commit.\n\n` +
    findings.map(f => `- ${f.file}:${f.line} [${f.severity}/${f.lens}] ${f.summary} -- ${f.failure_scenario}`).join('\n'),
    { label: `fix:round${round}`, phase: 'Fix', schema: REPORT })
  if (fix) { impl.commits.push(...fix.commits); impl.suite = fix.suite; if (fix.notes) leftovers.push(fix.notes) }
  if (round === 2) log('stopped after two review/fix rounds; remaining findings are reported as they stood')
}

phase('Check')
const check = await agent(
  `${WHERE}\n\nRun, in ${worktree}: \`${PY} -m pytest -q -p no:cacheprovider --cov=tanglepack --cov-report=json:.cand-cov.json\`, then \`${PY} tests/_tools/coverage_guard.py ${base_cov} .cand-cov.json --base-ref ${base}\`. ` +
  `Also report \`git diff --stat ${base}\` and the net added-minus-removed line count. Change nothing.`,
  {
    label: 'check', phase: 'Check', effort: 'low', schema: {
      type: 'object',
      properties: {
        suite: { type: 'string' },
        coverage_guard: { type: 'string', description: '"pass" or the newly missed lines' },
        diff_stat: { type: 'string' },
        net_lines: { type: 'integer' },
      },
      required: ['suite', 'coverage_guard', 'diff_stat', 'net_lines'],
    },
  })

return { commits: impl.commits, notes: impl.notes, leftovers, last_findings: findings, check }
