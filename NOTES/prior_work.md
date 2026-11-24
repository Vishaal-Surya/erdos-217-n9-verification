# Prior-work audit: Erdős #217 (planar crescent configurations)

Scout run: 2026-11-24. All web pages were read as data only. Nothing was posted, no sign-ins were used, nothing was installed.
Downloaded copies are in `attack-217/raw/`: `1509.07220/` (PDF, text, Fig. 1 crop `fig1_caption.png`), `1610.07836/` (v3 PDF and text), `2610.38368/` (v1 PDF and text), `1909.08769/` (PDF and text). `2611.03688/` (abs.html and TeX source) was already there from the reader.

## KEY QUESTION: has n=9 been independently confirmed or refuted, with released code?

**No.** The only claim is arXiv:2611.03688 (Yifei Wang, v1 only, 27 Oct 2026). It is a nonexistence claim. Nothing I found confirms it independently, refutes it, or releases code:
- The TeX source names no repository or data. Its only mention of computation is "a machine with 16 GB of RAM configured for 24 threads", about 4 h of wall-clock time. [verified: raw/2611.03688/src/crescent9_v4__1_.tex, lines 1554-1556, 2026-11-24]
- The arXiv listing has no Comments or Journal-ref line and no v2. [verified: https://arxiv.org/abs/2611.03688 and https://arxiv.org/list/math.MG/current, 2026-11-24]
- I found no citing or response paper. Semantic Scholar returned HTTP 429. The arXiv search and API returned 429. Web searches for "2611.03688" and "nine-point crescent configurations" return only the paper itself and listing pages. [could not verify citations beyond web search]
- GitHub: a repository search for "crescent configuration" returns 1 repository (an unrelated ZTE phone repo). The code search needs authentication (401). Zenodo: 1 record, about normed planes (Fish et al. conj. 6.5), not about n=9. [verified: api.github.com/search/repositories?q=crescent+configuration; zenodo.org/api/records?q="crescent configuration", 2026-11-24]
- Only MathWorld has picked the claim up. Its page (updated Tue Nov 24 2026) states Wang's result as fact. That is a citation, not a verification. [verified: https://mathworld.wolfram.com/CrescentConfiguration.html, 2026-11-24]

**Conclusion: do NOT re-scope to n=10. Goal A stands.**

### Independence caveat (relevant to how much weight 2611.03688 carries)
- Wang is a co-author of 1610.07836v3 (the classification paper). Chen (2610.38368) is at the same institution (Nankai, LPMC). Both papers thank Palsson. So the only "other" recent work comes from the same circle, and Chen's 3D paper (28 Oct 2026) does not cite Wang's n=9 paper at all. [verified: 1610.07836v3 author line; 2610.38368v1 p.1 and refs; Wang tex line 41]
- Wang's acknowledgment reads "OpenAI's GPT-6-Astra model assisted in the search for the constructions presented in this paper". That makes no sense in a paper that proves nonexistence. It is near-verbatim the sentence in Chen's 3D paper ("OpenAI's GPT-5.5 Pro assisted in the search for the constructions presented in this paper"). This suggests template or boilerplate copying, so read the paper carefully. [verified: Wang tex lines 1558-1563; 2610.38368v1.txt lines 330-331]
- 1610.07836v3 says "Computational source code is not included with this article" (Remark 1.4), and its certificates are "available from the authors upon request". It also says some figures and tables were generated with GPT-6. [verified: raw/1610.07836/1610.07836v3.txt lines 71, 159-160]

## Per-paper verification

| PROMPT.md claim | Status | Evidence |
|---|---|---|
| arXiv:2611.03688 exists: Wang, "no 9-point planar crescent config", 27 Oct 2026, single author, v1 only | **TRUE** | Title "The nonexistence of nine-point crescent configurations". Yifei Wang (Nankai). v1 Tue 27 Oct 2026 01:19:33 UTC, 29 KB. math.MG primary, math.CO cross-list. MSC 52C10; 05C30, 14Q10, 68R05. No comments field. [verified: https://arxiv.org/abs/2611.03688, 2026-11-24] |
| Wang: core = delete endpoints of mult-1 and mult-2 classes, giving a 6-point core with ≤6 distances | **TRUE** (as stated by Wang) | Tex lines 185-191: "Deleting one endpoint of each pair, and further points if necessary, leaves six points... It determines at most six distances." The Reader's job is to check the details (Gram minors, dimension 126, order ≤12, F/U). |
| Wang: "~4 h on 24 threads (garbled in the source)" | **Partly false**: the source is not garbled | "approximately four hours on a machine with 16 GB of RAM configured for 24 threads. This is the wall-clock time of that execution, not a general running-time bound." (tex lines 1554-1556) |
| Wang does not settle n≥10, because deleting a point does not preserve multiplicities | **TRUE** | Tex lines 179-182. |
| arXiv:1509.07220 (Burt, Goldstein, Manski, Miller, Palsson, Suh) | **TRUE** | v1 only, 24 Sep 2015, comment "Version 1.0, 4 pages". Published as Integers 16 (2016) A38 (per the bibs of Wang and Chen and per MathWorld; I did not open the Integers page). Theorem 1.3 gives crescent configurations in R^{n-2}. [verified: https://arxiv.org/abs/1509.07220, 2026-11-24] |
| Burt Remark 3.1: exhaustive search of a 91-point hexagonal patch, no n=9, "~900 CPU-hours" | **TRUE, with a wording fix** | The paper says "over 900 hours of computation" for the "naive implementation" on a "parallel computing cluster". It does not say "CPU-hours". [verified: raw/1509.07220/1509.07220v1.txt p.4] |
| Palásti n=8 coordinates as printed in Burt Fig. 1 | **TRUE, and the coordinates pass** | See the next section. |
| Constructions n=3…8: Pomerance n=5; Liu/Palásti n=7; Palásti n=8 | **TRUE** (secondary sources) | Burt abstract: Pomerance n=5, Palásti n=7,8. Wang (tex 104-110) and Chen (p.2): Liu and Palásti found n=7 independently; Palásti 1989 found n=8. I did not open the primary sources (Palásti 1987, Studia Sci. Math. Hungar. 22:447-448; Palásti 1989, Period. Math. Hungar. 20:231-235; Liu 1986, Math. Chronicle 15:29-33). [could not verify primary] |
| arXiv:1610.07836 v3, Oct 2026: n=4: 3 classes; n=5: 31 of 49 realizable; n=6: 6 of 6837 realizable | **TRUE** | v1 25 Oct 2016, v2 10 Jan 2019, v3 Mon 26 Oct 2026 20:53:55 UTC. Comment: "31 pages; updates made to algorithms and 6 points configurations added". The v3 authors are Chen, Durst, Huynh, Miller, Palsson, **Wang**. [verified: https://arxiv.org/abs/1610.07836, 2026-11-24] |
| The strictness filter can discard realizable classes, so 1610.07836 is not a complete classification | **TRUE** | 1610.07836v3 Remarks 3.5-3.6: Liu's parallelogram type is excluded because it is "not a strict distance class". Figure 6 says the non-rectangular parallelogram is not cyclic, so it is realizable in general position. Chen (2610.38368 p.2): "Strictness is an additional combinatorial restriction, not a necessary condition". Wang (tex 137-145) says the same. |
| arXiv:2610.38368 exists: Chen, 3D, crescent configurations in R^3 for every n via a helix | **TRUE** | "Crescent configurations in three dimensions", Yuda Chen. v1 Wed 28 Oct 2026, 8 KB. Theorem 1.1: P_j=(cos(j/n), sin(j/n), 2j/n), so D(n)≤3. It explicitly does not decide the plane. Acknowledges GPT-5.5 Pro. Note: the ID prefix is 2610 but the submission date is 28 Oct. That is consistent, and the ID is real. [verified: https://arxiv.org/abs/2610.38368, 2026-11-24] |
| teorth/erdosproblems lists #217 as open (synced 2026-11-23) | **TRUE** | data/problems.yaml: status open (last_update 2025-08-31), formalized yes (2026-10-23), forum posts: 3. The local clone HEAD cfd0fdb (2026-11-23 21:07) equals the live HEAD. [verified: raw/erdosproblems and https://api.github.com/repos/teorth/erdosproblems/commits, 2026-11-24] Oddity: the same entry has `formal_status: unformalized` but `formalized: yes`. |
| erdosproblems.com "froze statuses and proof claims on 2026-10-06" | **Partly verified** | The teorth README says the progress chart was frozen on Oct 6, 2026, under a policy announced at erdosproblems.com/forum/thread/blog:9, "under which that site no longer displays the statuses of problems, nor the count or percentage of problems solved". I found nothing about freezing "proof claims". erdosproblems.com returns 403 to automated fetches, so I could not read /217, its forum thread (3 posts per teorth) or blog:9. [could not verify the "proof claims" part or the live page] |
| formal-conjectures has an Erdős 217 file | **TRUE** (not listed in PROMPT; new) | `FormalConjectures/ErdosProblems/217.lean`, added 2026-10-23 in commit eace6d1 "Add Erdős problem 217 (#4064)" (only commit). It states `erdos_217.answer = {n | 0<n ∧ ∃ P, P.card=n ∧ InGeneralPosition P ∧ (distanceSet P).card=n-1 ∧ image multiplicities = Icc 1 (n-1)}` with `answer := sorry`, `@[category research open]`. It contains no n=9 lemma. [verified: local raw/formal-conjectures and https://api.github.com/repos/google-deepmind/formal-conjectures/commits?path=FormalConjectures/ErdosProblems/217.lean, 2026-11-24] |

### Other related work found (not in PROMPT)
- arXiv:1909.08769 (Fish, King, Miller, Palsson, Wahlenmayer, Crescent configurations in normed spaces, Integers 20 (2020) A96). It contains no 9-point Euclidean result, but notes that Palásti's 8-point configuration fits in a 37-point triangular-lattice patch. It has no n=7 coordinates.
- Clemen, Dumitrescu, Liu, Acta Math. Hungar. 177 (2025) 231-245, Obs. 5.1: arbitrarily large sets with these multiplicities that violate general position (cocircular). Cited in Chen p.2. [secondary]
- Zenodo 10.5281/zenodo.24384613 (2026-11-20): a counterexample about a normed-plane conjecture of Fish et al. Not relevant to n=9.
- primereboot.pages.dev/erdos/problem/217 labels #217 "SOLVED" with "last verified Jan 14, 2026" and gives no basis. It is an unreliable mirror, so ignore it.
- MathOverflow and erdosproblemaday.com: web searches found no page on crescent configurations or #217. Google Scholar was not queried directly (it blocks automated access). [could not verify]

## Known example coordinates

### n = 8 (Palásti 1989, as printed in Burt et al. Fig. 1)
The caption reads, from the rendered PDF (`raw/1509.07220/fig1_caption.png`): "Palasti came up with a crescent configuration for n = 8 [Pal89]. Points have cartesian coordinates"
(0,1), (√3,0), (2√3,0), (5√3/2, 5/2), (3√3/2, 9/2), (√3/2, 7/2), (3√3/2, 7/2), (√3, 2).
[verified: arXiv:1509.07220v1 p.2, 2026-11-24]. These match PROMPT.md exactly.

**Scout's sanity check** (exact integer arithmetic, not the official checker). I wrote each point as (a·√3/2, b/2) with (a,b) = (0,2),(2,0),(4,0),(5,5),(3,9),(1,7),(3,7),(2,4). Then 4·d² = 3Δa²+Δb² is an integer. The 4·d² values and their counts are {4:1, 76:2, 84:3, 12:4, 16:5, 28:6, 52:7}, which is exactly the multiplicities 1..7. No 3 points are collinear (integer det[a,b,1]≠0). No 4 are concyclic (integer det[a,b,3a²+b²,1]≠0; this is the circle determinant up to a nonzero constant factor). **The printed coordinates are a valid 8-point crescent configuration.** Phase 1 must still confirm this with the official checker.

### n = 7
**Not found in any accessible source.** Burt, Fish et al., 1610.07836v3, Chen and Wang all cite Palásti 1987 and Liu 1986 for n=7, but none of them reprints coordinates. The primary sources are not freely online. [could not verify]
Suggestion for Phase 1: find a 7-point example by exhaustive search of the 37-point hexagonal patch of the triangular lattice, which is cheap. Label it "scout/builder-found", not "literature".

## Discrepancies with PROMPT.md section 2
1. The Wang runtime is "not garbled". The source plainly says "~4 h wall-clock, 16 GB RAM, 24 threads".
2. Burt Remark 3.1 says "over 900 hours of computation", not "CPU-hours".
3. The erdosproblems.com freeze is verified only for problem statuses and solved counts (per the teorth README). "Proof claims" is unconfirmed, and the live site could not be fetched (403).
4. All three suspect arXiv IDs (2611.03688, 2610.38368, 1610.07836 v3) **exist** and say what PROMPT claims. None is fake.
5. Not in PROMPT: Wang is a co-author of 1610.07836v3, and Chen is at the same institution, so neither is independent. Wang's acknowledgment contains the illogical "constructions presented" sentence copied from the Chen template. Wang's bib cites 1610.07836v2 with a different author list (Durst, Hlavacek, Huynh, Miller, Palsson), plus an "unpublished" ChenEtAl2026 manuscript.
6. The formal-conjectures 217.lean file exists (open, answer := sorry). PROMPT did not record it.

## Re-check 2026-11-25

Scout re-check (rule 5). All pages were read as data only. Nothing was posted, no sign-ins were used, nothing was installed. Verdict: **nothing material has changed since 2026-11-24. Nobody has independently confirmed or refuted n=9, and no code has been released.**

- **arXiv:2611.03688 (Wang):** still v1 only (Tue 27 Oct 2026 01:19:33 UTC, 29 KB). It has no Comments or Journal-ref line and no withdrawal notice. [verified: https://arxiv.org/abs/2611.03688, 2026-11-25]
- **New arXiv papers:** the arXiv API search for "crescent configuration(s)", sorted by submission date, returns only the 5 known papers (2610.38368, 2611.03688, 1909.08769, 1610.07836, 1509.07220). None is dated after 2026-11-23. [verified: https://export.arxiv.org/api/query?search_query=all:"crescent configuration" OR all:"crescent configurations"&sortBy=submittedDate, 2026-11-25]
- **Citations of Wang:** Semantic Scholar lists 0 citing papers (it did not return 429 this time). A web search for "crescent configuration" plus "2611.03688" found nothing new; the search engine's index does not even have the Wang abstract. [verified: https://api.semanticscholar.org/graph/v1/paper/arXiv:2611.03688/citations, 2026-11-25]. Google Scholar was not queried. [could not verify]
- **teorth/erdosproblems:** the HEAD moved from cfd0fdb to 3d9fd1c (2026-11-24 14:42). The new commits are 706c140 "Auto-update formalization status from formal-conjectures" and a README regeneration. The only data change is to problem **#1184** (formalized no→yes), not #217. The #217 entry is unchanged: status open (last_update 2025-08-31), formalized yes (2026-10-23), forum posts 3, and still `formal_status: unformalized`. No commit message mentions 217 or crescent. [verified: https://api.github.com/repos/teorth/erdosproblems/commits and https://raw.githubusercontent.com/teorth/erdosproblems/3d9fd1c/data/problems.yaml lines 3954-3970 and 21688-21704, 2026-11-25]
- **formal-conjectures 217.lean:** it still has only commit eace6d1 (2026-10-23, #4064). `erdos_217.answer := sorry`, `@[category research open, ...]`, and there is no n=9 or other small-case lemma. [verified: https://api.github.com/repos/google-deepmind/formal-conjectures/commits?path=FormalConjectures/ErdosProblems/217.lean and raw main/217.lean, 2026-11-25]
- **GitHub:** a repository search for "crescent configuration" still returns 1 result (the unrelated ZTE Crescent Android repo). The code search was not attempted because it needs authentication. [verified: https://api.github.com/search/repositories?q=crescent+configuration, 2026-11-25]
- **Zenodo:** still 1 record (10.5281/zenodo.24384613, about normed planes, 2026-11-20). It is not about n=9. [verified: https://zenodo.org/api/records?q="crescent configuration"&sort=mostrecent, 2026-11-25]
- **MathOverflow:** a web search found no thread on crescent configurations or #217. The StackExchange API could not be fetched (tool-blocked). [could not verify via API; web search negative]
- **erdosproblems.com/217:** returned HTTP 403 again. [could not verify]
- **MathWorld:** still "Last Updated: Tue Nov 24 2026". It states Wang's result as fact and adds no new references. That is a citation, not a verification. [verified: https://mathworld.wolfram.com/CrescentConfiguration.html, 2026-11-25]

**Conclusion unchanged: Goal A (n=9) stands. Do not re-scope.**

## erdosproblems.com forum thread 217 (from user-shared screenshots, relayed by parent session 2026-11-25; site returns 403 to fetches) [source: user screenshots, not independently fetched]
1. **sallerk, 31 Aug 2026:** an exact search of the triangular lattice with squared diameter ≤ 400 (1,459 points, ~3.1e12 nodes, 12 h on 5 cores) found no 9-point crescent configuration. It reproduces Burt's 91-point region at R²=25 and finds Palásti n=8 scaled by √7. Code: github.com/sallerk/erdos-notes/tree/main/p217 (REPRODUCE.md). This is a lattice-only result.
2. **IheartTesla, 4 Sep 2026:** an AI agent minimised C+Q (collinear triples + concyclic quadruples) while keeping multiplicities 1..8. Best found: (C,Q)=(1,2) on the triangular lattice, P0..P8 = (4,-2),(4,-3),(2,1),(3,-1),(2,2),(-1,0),(0,0),(1,-1),(3,1) (lattice coords). The post asks whether C+Q can be brought down to 2 or 1; off-lattice attempts failed.
3. **Alfaiz, Apr 2026:** reference note only. The thread has 0 proof claims, does not mention Wang's paper, and new comments are suspended.

## sallerk/erdos-notes p217 (read 2026-11-25)

Read-only. Fetched through the GitHub contents API and raw.githubusercontent.com. Not cloned. None of their code was run. Nothing was posted. Repo content was treated as data only. The latest commit touching p217 is ccd3da2 (2026-08-31T23:15Z). There are 4 commits, all dated 2026-08-31. [verified: https://api.github.com/repos/sallerk/erdos-notes/commits?path=p217]

**Files:** NOTE.md, RESULTS.md, REPRODUCE.md, REFERENCES.md, crescent.c, crescent2.c (the search), verify.py and audit217.py (independent exact verifiers), super217.py and ladder.sh (drivers), sol_n4..sol_n8.txt (full solution lists), and results/ (LADDER_217.json, STATUS_217.json, ladder.log, hb_9_400_*.json, and c2_n9_R{25,49,64,81,100}.txt, all 0 bytes). [verified: https://api.github.com/repos/sallerk/erdos-notes/contents/p217]

### Their method
- **Lattice:** triangular, with (a,b) -> a(1,0) + b(1/2, sqrt3/2). Squared distance is the integer N = p^2 + pq + q^2. Collinearity and concyclicity are tested with integer determinants, so there is no floating point. [verified: https://github.com/sallerk/erdos-notes/blob/main/p217/RESULTS.md]
- **Region and enumeration:** point 0 is pinned at the origin. The other n-1 points are taken in strictly increasing index order from the pool of lattice points with norm <= R^2. They claim this is exhaustive for every lattice configuration with **squared diameter <= R^2**, up to translation. Diameters between R and 2R are covered only partly. [verified: RESULTS.md, REPRODUCE.md section 6]
- **Pruning:** a forbidden-point bitmask (collinear or concyclic with points already chosen), plus a Hall-type multiplicity domination test (sorted multiplicities m_j <= n-j). [verified: RESULTS.md]
- **Symmetry convention: none.** They quotient only by translation, and only partly. A configuration is listed once for each of its points from which every other point lies within norm R^2, so the counts are raw and are **not** reduced by the 12 lattice symmetries. [verified: RESULTS.md; consistent with the sol_n*.txt contents, where every line starts at (0,0)]
- **Multiplicities:** they use the multiset {1..n-1} and do not require them to be monotone in the distance. This is the same as our definition. [verified: https://github.com/sallerk/erdos-notes/blob/main/p217/NOTE.md]

### Their numbers
- **Positive controls at R^2 = 49 (187 points):** n=4: 2892, n=5: 4116, n=6: 8904, n=7: 2964, n=8: 156 solutions. Times were 0.11 / 0.32 / 1.53 / 5.78 / 20.54 s. The n=8 run took 161,595,043 nodes and 15.3 s in REPRODUCE.md. The line counts of sol_n4..8.txt match these numbers exactly, and every point in those files has norm <= 49. [verified: RESULTS.md, REPRODUCE.md, raw sol_n*.txt]
- **n=9 ladder, 0 solutions on every rung:** R^2 = 25/49/64/81/100/144/196/256/324/400. The pools have 91/187/241/301/367/517/721/931/1165/1459 points. Nodes run from 16,636,430 (R^2=25, 2.8 s) to 1,826,221,345,879 (R^2=400, 25,082 s). The total is 43,203 s (12.0 h) on 5 cores and about 3.07e12 nodes. LADDER_217.json and STATUS_217.json both say COMPLETED with 5/5 shards for rungs 100-400 and HEALTHY=true. Rungs 25-81 appear only in the markdown and in ladder.log, which was cut off mid-run at R^2=100. [verified: https://raw.githubusercontent.com/sallerk/erdos-notes/main/p217/results/LADDER_217.json, STATUS_217.json, ladder.log]
- **Burt et al.'s 91-point region:** they report it is identical to the R^2=25 pool. Their 0 at n=9 on that region reproduces Burt Remark 3.1 in about 3 s. [verified that they state this: RESULTS.md. Their set-equality check lives in audit217.py check 6, which was not run → could not verify]
- **Palasti n=8:** they say their first n=8 solution is Palasti's set scaled by 7 in squared distance, with an explicit bijection (7,4,6,5,3,0,2,1). [verified that it is stated: commit 91fd273 message and REPRODUCE.md. audit217.py was not run → could not verify]
- **Reproducibility:** they released the full source, build commands, the outputs each command printed, the complete n=4..8 solution lists, and per-shard records. In principle that is enough to reproduce, but we have not done so. Their ladder.sh originally depended on a local env.sh that was missing from the repo; it was fixed in ccd3da2. [verified: REPRODUCE.md and the commit message]

### Comparison with ours (37-point radius-3 hexagon, up to 12 symmetries + translation)
Their raw counts (2892, 4116, ...) are **not directly comparable**. They use a different region (a norm disc pinned at the origin, not a fixed patch) and no symmetry reduction.

They **are** comparable once reduced, because their files contain everything we need. Any configuration that fits in our patch has squared diameter <= 36, which is at most 49. So it appears in their R^2=49 lists: translate any of its points to the origin and every other point has norm <= 36. I reduced their published sol_n*.txt files with my own script, which shares no code with theirs: canonical form under the 12 symmetries (rotation (a,b)->(-b,a+b), reflection (a,b)->(b,a)) plus translation, keeping only orbits that some translate fits inside max(|a|,|b|,|a+b|) <= 3. Each kept representative was re-checked in exact arithmetic (no 3 collinear, no 4 concyclic, multiplicities = {1..n-1}). None failed. (The scout's scratch script was not retained in this release.) [verified: computed 2026-11-25 from https://raw.githubusercontent.com/sallerk/erdos-notes/main/p217/sol_n{4..8}.txt]

| k | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| ours (37-patch, orbits) | 24 | 36 | 19 | 26 | 8 | 1 | 0 |
| theirs, reduced to the 37-patch orbits | n/a (they do not search n=3) | **36** | **19** | **26** | **8** | **1** | **0** (implied) |
| theirs, all orbits within R^2=49 lists | n/a | 187 | 123 | 231 | 72 | 5 | 0 |

- **k=4..8: exact agreement.** These are two independent searches.
- **k=8:** the single in-patch orbit from their data is ((0,0),(0,2),(1,-2),(2,0),(3,1),(4,-3),(4,-1),(5,-1)). By count it should be Palasti's, matching our known_examples.json. I did not check point-for-point that it is our stored representative. [could not verify identity, only the count]
- **k=9:** our patch's squared diameter is <= 36, inside their exhaustive R^2 <= 400 claim (R^2=49 rung already suffices), so their 0 implies our 0. They agree. [their claim verified as stated; their search not re-run]
- **k=3:** they do not search n=3, so there is no comparison. [could not verify]
- The "all orbits" row only counts orbits that contain a point seeing all the others within norm 49. It is not an exhaustive count for squared diameter in (49, 196]. It is shown for scale and should not be quoted as a count.

**Bottom line:** this is independent confirmation of our k=4..9 patch counts. Their n=9 negative result (lattice, squared diameter <= 400) strictly contains our patch search, so our k=9 = 0 is not new information. It is a consistency check.
