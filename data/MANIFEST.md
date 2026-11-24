
## Phase 2 pattern enumeration (Builder E, 2026-11-24)

Producers: src/enum/PatternEnum.java (mode A), src/affine.py, tests/soundness_run.py. See NOTES/phase2_counts.md.

| file | lines | sha256 (file as stored) | sha256 (uncompressed) |
|---|---|---|---|
| patterns3_local.txt | 3 | 7e087d0e8410b888806587482cb782cd221fade886a310acedda5d60d2865c5c | (same) |
| patterns4_local.txt | 22 | 3b31c7ed170f32fb5738abd3538db1a7f6fb0036c46f3636f6f21081d8d926a0 | (same) |
| patterns5_local.txt | 932 | 09c4602042ae14a962f9c8173b4f6bd3f99ca7b10085896a0293ac94d3ba531e | (same) |
| patterns6_local.txt.gz | 625472 | 7bc1ef99164999dacb385bff9fbba9e4c90d402805aaaeb1aeecefe3f2e3c3fe | 1df5ffc01c74da3ed53c1563b7b3a886278f330dc77a2e93bccc7ef2b0bc16f6 |
| patterns6_affine.txt.gz | 344835 | f15fd555692c94cafee3fd2251f6849afa70130cde59a277881e33558f643869 | ae330ae246b12f2f12515d6093fc5d0756dd3f1b6a9da1f500a6119a4fac9bc7 |
| soundness_realizable6_W8.txt.gz | 202 | dcd966bc75cc61adbfd0a861c9b227534aab815db978193bb0b9a62ad9e82bfa | 4ba9db58beeb35d5d8d31a1b351b8d1c8e30f71df78a8c867ea5e29828215a9a |
| soundness_realizable6_W10_111_101_301_103.txt.gz | 194 | febe42764aab5ca2b92d4aaacf1888a029f6987166105cb5c424ae6d1b7f153a | 7082829b265905be775658bb91709f2e19d9184cd2cdf3c4b91d8f46ebc996e9 |

## E1 full run (orchestrator extract from data/algebra_E1_chunks/, 2026-11-24)
- data/E1_F.txt.gz  sha256 6a6b6e1342cc54a527eda95609c04f9fd5ffcbe9525e631cce1b29ef6095f9fe  (4,450 lines)
- data/E1_U.txt.gz  sha256 cbd6781f94f2e2f9de47202ee7b60965e7c5f69e977789accee356f376e4214a  (187 lines)
- data/E1_excluded_certs.jsonl.gz  sha256 7f57361dda5533a2fd55696df84323c02d6a735e4c8344c1ea72c0797a54e8fe  (340,198 certs)

## Builder D: F reconstruction + extension (2026-11-24; see NOTES/extension.md)
| file | sha256 (as stored) | sha256 (uncompressed) |
|---|---|---|
| F_recon.jsonl.gz | e54f6dd29f93816bd6317df3228798f63142dae69724710495c9c9e655e5de88 | 685c7d3f10bb2fe0eb79f934ca3c32b90d62ad61df5128f0efba15715edf8f28 |
| F_recon_redo.jsonl.gz | a77e099d87f4e992485f56f93a1ed542ac38b265b93d1881602488728b98de15 | 99a55f7512cf9c3c5e396b6cdcc9fb76685bdbf5801816fa2e288d69050c9df3 |
| F_recon_final.jsonl.gz | 67604ae28499a1f548635a662394938519ddd1b1cfd86945a9e4543139e29a29 | ccb936cd242bdb07c1a6b60121ab03510719c72d07d21085d9e59afb46a3fa6c |
| F_extend.jsonl.gz | 0c38488f445362256d35658303a24ae7bce19cb9822e55d4e25b8060c1d1079b | a06bf92b62f3f0fd9677831b54f75955a8aba32de062e4298e1e42e3b252c22b |
| F_extend_full.jsonl.gz | 7e8db5a608ecd89833803828904f06ca96bbb9e37ec47a63a4796e0b175c6194 | 442e13d40f7e8241278d29082a95f2202a7d915e2933865a639b3a5a0c40a8bb |
| F_recon_e1fallback.jsonl | 6efd3ad05e246c19f88df5f8fadfcf4a12bb5967888945bdd67c752341c4e99c | (same) |
| F_recon_redo_in.txt | 3ad1ea32b2ac72a6831e6db257cb4534062083cd6eae84acbbf3b1c00f2a5519 | (same) |

## Coverage search (Builder C, 2026-11-24)

Producers: src/coverage/CoverageSearch.java, affine_n.py, verify_survivors.py, run_real.sh. See NOTES/coverage.md §5. Input U = data/E1_U.txt.gz.

| file | lines | sha256 |
|---|---|---|
| E1_U.txt.gz | 187 | cbd6781f94f2e2f9de47202ee7b60965e7c5f69e977789accee356f376e4214a |
| coverage_E1_run.log | 20 | 5c289b85a17b6e28b8704aa1aaffce3d1ab989e0691fa1ba6c4627fbfe9134ef |
| coverage_E1_found.txt | 487 | 5b5b4cd68a9ccb599d2cd4af9f9164db2b7f8cbb0ad3fff23767b0cbdbdfa644 |
| coverage_E1_affine.txt | 487 | 2ebd2c8c18756b5b45e9805a8f6aa44071d4db0ce63bc82c48c7a3589cc06dc4 |
| coverage_E1_level7.txt | 13 | aa9093c7cf3b21d0ce5fe96465bf7ebeef7dc78e249aa94f55a48a453f8ae814 |
| coverage_E1_level7.log | 5 | b3266158e86f92f8901c93d3c7e0ae9bd369af9fac144ae0a5af26aca9a5dc55 |
| coverage_E1_run_paranoid_spot.log | 4 | 8e9c0c138b999ccba50a877258cc1e8bd7cbef5a4170e222080c85e1203c073b |
| coverage_planted_s1.log | 17 | b908175fc63fc7bfc861f063cc649e5740891bf571eb5e3f98c75df7c614b2e9 |
| coverage_planted_s2_noise300.log | 12 | b584bc1c941cf56bc053b881ad1ed3570a2ed67029d81334197746fc69dccaa2 |
| coverage_standin_affine300.log | 10 | 16c370aa69a26deefbcf59543bdac8b9b4418d139670932112b2f8b716324f7f |
| coverage_n7_localU.log | 26 | 6d64628d49a1461e201d779bf4cc11cc1dee39e40104ced157f804ad6488576e |
