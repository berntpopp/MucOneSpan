# Consensus anchor configuration fixes

2026-09-14. No final seeds or calling were used.

## Contract and causes

Anchor components now contain up to the configured number of available bases, and coordinates subtract the actual sliced length. Oversized anchor requests therefore cannot create negative coordinates through subtraction of unavailable bases. The private `_find_anchor` API still expects an anchor start; its callers translate the expected VNTR boundary into the corresponding start. Tolerance is now the inclusive permitted boundary displacement, including zero for an exact boundary. Left-flank anchors use the suffix of the prefix actually selected by the ladder, and right-flank anchors use the corresponding selected prefix. Ladder construction itself was not modified.

The corrected left anchor legitimately changes default output for some existing consensuses: a one-base insertion in the actual 500-base left flank previously left an extra flank base in the VNTR. This is a pre-freeze correctness correction, not a threshold adjustment or evidence of improved biological recovery. The previous default parity claim does not apply universally after this correction. Correcting the search center also changes the accepted displacement interval from an anchor-start-centered interval to the documented boundary-centered interval.

## Red to green and checks

Before implementation: five newly exposed failures, five passing controls. Failures included anchor size 20 producing trim_start=-12 despite an accidentally correct Python negative-index slice, anchor size 100 producing an empty VNTR, zero-tolerance exact anchors marked absent, and short/default flank prefixes retaining an extra base after insertion. After implementation and additional symmetric displacement controls: 38 unit tests passed with no skips. Twenty real bcftools/samtools integration tests passed with no skips. Ruff and mypy passed. The specifically authorized legacy coordinate fixture in test_consensus.py now builds its reference with the same left prefix as the ladder.

```sh
uv run --locked --all-extras pytest -q tests/unit/test_consensus.py tests/unit/test_consensus_settings.py --no-cov
PATH=/home/bernt-popp/miniforge3/envs/env_clair3/bin:$PATH uv run --locked --all-extras pytest -q tests/integration/test_phase_consensus.py tests/integration/test_variant_concordance.py --no-cov -rs
uv run --locked --all-extras ruff check src/muc_one_span/consensus.py tests/unit/test_consensus_settings.py tests/unit/test_consensus.py
uv run --locked --all-extras mypy src/muc_one_span/consensus.py
```

## Existing cached full-consensus comparison

Selection is the explicit 77-input development inventory in `phasing_debug/coverage5/inventory.json`; paths were mapped to the original `development/{group}/{sample}` coverage-10 results, not the coverage-5 reruns. Both expected allele slots were enumerated even when absent. Each available full consensus was trimmed with the saved preconfiguration source and the corrected source using the same bundled dictionary, flank length 500, and default settings. Before comparison, both actual per-allele reference flank prefixes were checked against the dictionary. Existing result files were read only; temporary output files were deleted. SHA-256 values below are of raw DNA sequence bytes for trimmed sequences and file bytes for full-consensus/reference/source files. All changed sequences are listed; unavailable slots remain explicitly accounted.

```json
{
  "inventory_sha256": "66182019aba7b89c9d0d6317be714fe7e25585f7f1c2eff8d94f3d0b007364ac",
  "old_source_sha256": "d6cca7ca2a4ec7070069b1d0c2b2448947f9790db8c05b8b33d93f6fc6022e0d",
  "new_source_sha256": "9a54eed7ba4be57ef42d828b973f2cf855cd0fb6c8e35eeb1a084555eccdc305",
  "expected_inputs": 77,
  "expected_allele_slots": 154,
  "missing_slots": [
    {
      "group": "heldout",
      "sample": "sample_heldout_normal_60_62",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n10_seed1701",
      "allele": "allele_1"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n10_seed1701",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n10_seed1702",
      "allele": "allele_1"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n10_seed1702",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n10_seed1703",
      "allele": "allele_1"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n10_seed1703",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n20_seed1701",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_dupc_60_80__n20_seed1702",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n10_seed1701",
      "allele": "allele_1"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n10_seed1701",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n10_seed1702",
      "allele": "allele_1"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n10_seed1702",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n10_seed1703",
      "allele": "allele_1"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n10_seed1703",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n20_seed1701",
      "allele": "allele_2"
    },
    {
      "group": "perturbations",
      "sample": "sample_normal_60_80__n20_seed1703",
      "allele": "allele_2"
    }
  ],
  "counts": {
    "available": 137,
    "changed_sequences": 0,
    "changed_contexts": 137,
    "old_matches_cached_trim": 137
  },
  "groups": {
    "hifi": {
      "available": 88,
      "changed_sequences": 0,
      "changed_contexts": 88
    },
    "ont": {
      "available": 6,
      "changed_sequences": 0,
      "changed_contexts": 6
    },
    "heldout": {
      "available": 11,
      "changed_sequences": 0,
      "changed_contexts": 11
    },
    "perturbations": {
      "available": 32,
      "changed_sequences": 0,
      "changed_contexts": 32
    }
  },
  "all_reference_flanks_verified_as_actual_ladder_prefixes": true,
  "rows_sha256": "411c46aa5494e047bf4b598e95410a718e4ea648762eb183e39f5425b99a1073",
  "changed_rows": []
}
```

## Per-sequence hashes and metadata difference

All 137 available cached outputs had identical old/new FASTA bytes and identical trim_start/trim_end coordinates. The only context change was `left_trim_method: fixed_anchor_not_found → exact_anchor` in all 137. The 17 unavailable slots are listed above and were not assigned parity results.

SHA-256 of the ordered per-sequence identity/hash table as canonical JSON (`sort_keys=True`): `fc6ee556505f72c49069d9257baec9024182cc4c5164149289cde02bb0476b4a`.

| Development input / allele | SHA-256 shared by old and new trimmed sequence |
|---|---|
| hifi/sample_asymmetric_25_140/allele_1 | `6d5e003f4697ab9e0335ca695990e9ce104549a45f760233056653583669c2b6` |
| hifi/sample_asymmetric_25_140/allele_2 | `7b84233007c700490393f0b37f7d3177863b202515bc98fc9f87a76f318b87d5` |
| hifi/sample_bench_5000/allele_1 | `5f92fb25ea4ae7a82d8cf2af41bcf9e0552acf486608b325c4bc65a371369076` |
| hifi/sample_bench_5000/allele_2 | `a83e8f8a1c4bd1986d2d8df86ab4cf19b309cfadcc5aa058ab7838dc88995726` |
| hifi/sample_bench_5001/allele_1 | `5109566ceb019347eb4397ec2850654ced9f6b57a44a6e4e0e9c4f46dd45ee7a` |
| hifi/sample_bench_5001/allele_2 | `aa47b50d41a4b872bbe623b799d28a9583db064891b1f680cc7a1ca9adaa513d` |
| hifi/sample_bench_5002/allele_1 | `a92d89767ead23fd05fc9e90c3e0c90cdc18ca324994b445c123858d7c2fc718` |
| hifi/sample_bench_5002/allele_2 | `fc4f376dbc0e3d0719a580c9de2a2adfe85664ebb5df8c723022333aa3647266` |
| hifi/sample_bench_5003/allele_1 | `6db6f1a30c62bceb9c5c3890225a4e78d22af0c92d32e0ee2613ace5f3fc9bd8` |
| hifi/sample_bench_5003/allele_2 | `8b28af1e762a5e01d220561b8f78bacd5f22d141f34bbe1ea715d4794bafa400` |
| hifi/sample_bench_5004/allele_1 | `14eb5ae5d27f23220fdfab07c565940525156b40b4c5f1d160a8a7969943b641` |
| hifi/sample_bench_5004/allele_2 | `818e28ed05e0f82d54e8f98b824fbfad176395b62e7641f0d6dd2684b9760781` |
| hifi/sample_close_51_58/allele_1 | `a92d89767ead23fd05fc9e90c3e0c90cdc18ca324994b445c123858d7c2fc718` |
| hifi/sample_close_51_58/allele_2 | `fc4f376dbc0e3d0719a580c9de2a2adfe85664ebb5df8c723022333aa3647266` |
| hifi/sample_del_60_80/allele_1 | `270d4a753e88dab0a64dd4c4c84c89af6399af83bdd7bc3a0c4d735cb80ef3af` |
| hifi/sample_del_60_80/allele_2 | `d5d0f03f91b980e0380c43a7fec50eec74fa8440ecb092f34591f724017cc1cf` |
| hifi/sample_dupa_100_120/allele_1 | `d9e47d4884e4302075d50bf1d15e3bb82667fefb577ed523335a3463393ea1ca` |
| hifi/sample_dupa_100_120/allele_2 | `6f190a96dbdc1108dccd6f15bc404f11f6896802ca8b5ed2efacada6ab9f2605` |
| hifi/sample_dupa_60_80/allele_1 | `1431eb73c00e91f6b0ff893674daabc1b0090695c160d222454e0793bb898907` |
| hifi/sample_dupa_60_80/allele_2 | `77ba07b9fb8de4468d7f1d0d1700481d5d1575be5d3adff7fdc7d4589374d5ca` |
| hifi/sample_dupc_100_120/allele_1 | `3e6bd31e0d0433920e962c525ae5add366502dd2462c1629750716b3cdbffbdf` |
| hifi/sample_dupc_100_120/allele_2 | `c00f49bcd0ded877a4474be20aadaf7da039cdf6c865f1146d2950a6e5c9f1ea` |
| hifi/sample_dupc_40_50/allele_1 | `36aa965cc51bdbf8039ffb7d610ea7f913b5895aa43881adf276533aa346d530` |
| hifi/sample_dupc_40_50/allele_2 | `c07d413b165f6e43e7f444a1568fd2642676b624ef0220a2916e7a4385155dbe` |
| hifi/sample_dupc_50_55/allele_1 | `cb40c8cbc2b21d9ed5682efee1bd0c832e080ec6cddcd6447d041fa32234d4a2` |
| hifi/sample_dupc_50_55/allele_2 | `e9b8abd0ad91dbcf1b2cd958dca9d0d4193050d9669716e156fca9a8ffdc0b81` |
| hifi/sample_dupc_50_57/allele_1 | `aa1f2e1fbeede0ecb374540ae183a3e801cc8e55071db8216ed2c7cf803cd5c9` |
| hifi/sample_dupc_50_57/allele_2 | `ba55a088eddde30c891506a5b2ac01c5dd99a6f2413d83acb0bbd74088775634` |
| hifi/sample_dupc_50_60/allele_1 | `ba2628790ac690bfb267e74f94decc2c2f39128067c26a3470e9dce4f679778d` |
| hifi/sample_dupc_50_60/allele_2 | `f22bf81899e94b7bda83d3e883ee6e2db39493cbd0f09e1e4c29ce52bdb4ae0b` |
| hifi/sample_dupc_60_80/allele_1 | `1389fb8cb0b5001e2226f6d47845682978711cc78c4c37b6ff568b6d258fb3a7` |
| hifi/sample_dupc_60_80/allele_2 | `f10fed38d9604b7acddef646160151b952e738ac6039375c7fbb53035fb380c8` |
| hifi/sample_dupc_60_80_cov50/allele_1 | `10b52bea0b5e173ea78e94a4763c2d2fae08820d5e47fb5abddd8bfbfda20e7a` |
| hifi/sample_dupc_60_80_cov50/allele_2 | `119086e323a17f6b508be7cdb40c533629633b1f1c199c250077c093820ea228` |
| hifi/sample_dupc_60_80_s2/allele_1 | `cd45f20bd9781b8524c5932767d1c7223f45d67cfed51a9d2c730f1abc954424` |
| hifi/sample_dupc_60_80_s2/allele_2 | `e6c53c1a667e631f6be0085e06ecd394337afebaf89d11a6cbaee3163819c7af` |
| hifi/sample_dupc_60_80_s3/allele_1 | `673fc1f7a3fc5c6377f2a508ed4ae47c5e7ee627b2394ddc36d4390b8053bd14` |
| hifi/sample_dupc_60_80_s3/allele_2 | `ee59b18c0e9d8efff1abba436614c3aafd63c64dda9712568ccccdbecaff1724` |
| hifi/sample_dupc_60_80_s4/allele_1 | `ec0d33a2c03ab7f89cc948761c31ad8b7df57755686d83af9cc2cf9b64f67f80` |
| hifi/sample_dupc_60_80_s4/allele_2 | `e4fc0cacd4bf779d02c5d02bdb4cbbbd6feca3a6607285e45bfe0c07b910efd0` |
| hifi/sample_dupc_60_80_s5/allele_1 | `7d704ca6307e616a8a3be69d39f5c146c9759622708a5acdfed275a473281d38` |
| hifi/sample_dupc_60_80_s5/allele_2 | `9fed1b4fb33f82115f3ad966cfce6ce41b8ac3f1c93b5097f813f5698ea72517` |
| hifi/sample_dupc_80_100/allele_1 | `65a607b80c9e0fde5290383b3d1be35d4c037412b9db73448cb31e9944960cf5` |
| hifi/sample_dupc_80_100/allele_2 | `d4196f47805d7f5d90072de0b0599b39bf91254126d18d719e0d667d5881e8cd` |
| hifi/sample_dupcccc_60_80/allele_1 | `89a1dc74994a98e5f38f623a3a0b628aee9b0933c62e1ab6567aad2e456870e7` |
| hifi/sample_dupcccc_60_80/allele_2 | `663b30e9457de425a5eb88bd8183aa1a2b8d878289a9d34b7ac0c2ef81f6a323` |
| hifi/sample_gap10_40_50/allele_1 | `e56174ee2e1d89697d7f1f6db5a267464a3d89a2a502fba21906816d2c25d7c3` |
| hifi/sample_gap10_40_50/allele_2 | `8e2421ae630766cbcf5c38bd00917caf83903a505b1e1b7986765fefdb3b685d` |
| hifi/sample_gap12_55_67/allele_1 | `2ff82604b7c04027dd35016d56238650462200ad0f91ae8203c2535723f906a1` |
| hifi/sample_gap12_55_67/allele_2 | `44910fcee1b6c50db94964afd5c49eec7ff1738cbc0322739bc4fc61b46cb9e7` |
| hifi/sample_gap3_50_53/allele_1 | `0fb1501c56a2ca4adf2ef274af51cd3e25e35c8a64dc73546b1702cf0f13ef24` |
| hifi/sample_gap3_50_53/allele_2 | `febb4ccb6f4091aa9342b2c176474e71c4b4be482159f9ad437c4d2b4418184b` |
| hifi/sample_gap3_80_83/allele_1 | `c8a286182ffd7b789ea7a47e4e57c11c751dff51b94380cf9eec45b5778c7676` |
| hifi/sample_gap3_80_83/allele_2 | `f62170acf78120cfa40d35011affacaed39257ff121ab0509b268174598236ae` |
| hifi/sample_gap4_40_44/allele_1 | `dbb89e0a3299474674101764ae8c8ae2d5e48d1d48b3dba0c7655f0e88e93e39` |
| hifi/sample_gap4_40_44/allele_2 | `63f508438770ae4fc030e43a460fa8890f8a5ef81632641b7cf463e65e58b5c8` |
| hifi/sample_gap5_55_60/allele_1 | `d21bdc7c2282cc4852da4861d7868de8eceb66f784d09eba6b18e6f67bbb6a57` |
| hifi/sample_gap5_55_60/allele_2 | `e893e7a74619409b14305ef9e71f1f26c0c4037656085a7123a23df39d04f71d` |
| hifi/sample_gap5_70_75/allele_1 | `c0123fbe035d1fa102dd16ee9793a6c94d9e578c78f4633af7f3d34815aac0af` |
| hifi/sample_gap5_70_75/allele_2 | `b39769f50ed5888c97d66acb08bdfb844ddb0e9f760427ea6cfdd0bf0f6f1071` |
| hifi/sample_gap6_30_36/allele_1 | `bac7f702e7944101acb670266466b1d5ff9bf100715df939b75e4ac16b139acb` |
| hifi/sample_gap6_30_36/allele_2 | `601342cee4209361a32786fdd95d56874ec6b27efbb8f020b9619fde5243abd2` |
| hifi/sample_gap6_45_51/allele_1 | `e93f1fb47d7b3f2ee16e981d63ea9a54390522399844454a02988ca75002e187` |
| hifi/sample_gap6_45_51/allele_2 | `fb227cb4c20468991c0264e6e9c2c078bcf40da8d5b69e98d5183cf07e05ca66` |
| hifi/sample_gap7_60_67/allele_1 | `9d909f59375501fadbace2eb9266d572080da4ff77d0f64e7cbfdaa4868c3781` |
| hifi/sample_gap7_60_67/allele_2 | `3fc0fd1f95c7444b3306c207c5c8765d85ddb8fca0f13e1b830cd7ea92cf6384` |
| hifi/sample_gap8_35_43/allele_1 | `cacd5f84cf58ecf52741956c1142abf84341ea2a4aa389ce5302985a09d26943` |
| hifi/sample_gap8_35_43/allele_2 | `b7fb13f0a7bd9f41eda942c9bea4275ed8c60c380f95be52abb00b980a5d9f5a` |
| hifi/sample_gap9_50_59/allele_1 | `331c5b3213651b4d1b4c618e78074fa4255d0c6762fa3ceee058e76998570ad8` |
| hifi/sample_gap9_50_59/allele_2 | `4ed0525da0b9f1696188c8bb2d126932dbe3a38a8d5b16100bdf07652f107ede` |
| hifi/sample_homozygous_60_60/allele_1 | `0bc9f181b1887836af8892ec35003f67531d5113d4b86710be9697717395d93e` |
| hifi/sample_homozygous_60_60/allele_2 | `87d23c2402d5d37ca4eac2131b61294e1cffcc0ec5b58a57252f77965746231e` |
| hifi/sample_insg_100_120/allele_1 | `1fc6adc78d6c1699bc47043fc333d0799067857124547d799609f13ec026c1d2` |
| hifi/sample_insg_100_120/allele_2 | `1884cf67789d77176b03bcead67369abbdf733e931fb8fdcd56003e17f513c98` |
| hifi/sample_insg_60_80/allele_1 | `d129953b339a94d929f9604fa95f3c163565617d0ecfdb4104c3d8c26f226859` |
| hifi/sample_insg_60_80/allele_2 | `7e2a7506c885da9eb7aa8d53dc21d9d255fae93a05fe921bcea5f26c79c09460` |
| hifi/sample_long_120_140/allele_1 | `12d7cc2db45cde1b6a6afa07e9dfb574ab8352049670e3d2074ba2752176ccd5` |
| hifi/sample_long_120_140/allele_2 | `20e66bc4378b32362c0cdd36501b89ea2ca9b6d181c12950644523b4fae3bd37` |
| hifi/sample_normal_100_120/allele_1 | `aed26998c2771c104160a2dda55c8f056461e5152c70bbd76b5ebf2e1650a429` |
| hifi/sample_normal_100_120/allele_2 | `42e6ad6b6cc3736c02f5ef1355262d63a6706de4973f0330302d788c10f3c1d3` |
| hifi/sample_normal_50_55/allele_1 | `a937de737c183a4aed678de534b86f0328e3851318df652fbfa723a1f04666a5` |
| hifi/sample_normal_50_55/allele_2 | `00f5a3bf8478066346bab29dbd157c12d66a95a232e0c92e786a47cded55ce98` |
| hifi/sample_normal_50_60/allele_1 | `fd8f0aebb0a78c78add5b8ac8c8a76843c43ea193dfbe54cc34775d247e13d4b` |
| hifi/sample_normal_50_60/allele_2 | `25e9fdc28b2a483985ea3e21d6722756dd206f3e0f224bb8a6f5ed35c4c40c50` |
| hifi/sample_normal_60_80/allele_1 | `372fb6f499ed4f4f07bb6274e5b2e8c1060bb6fcefec8311aecdef3b6d3661c8` |
| hifi/sample_normal_60_80/allele_2 | `86faafad8d711c69b1664314ca88f29e0414bf2cd88e003322b3ab767fc479ca` |
| hifi/sample_short_25_30/allele_1 | `81aba661fbd8482cb0db811bb2b9ac6ad6812697b42d62051f520464f97a8de6` |
| hifi/sample_short_25_30/allele_2 | `6a41e653fa9a620f589f6d2b4045f66e03b17b847141a84f64dbf5bf534099a4` |
| ont/sample_ont_dupa_60_80/allele_1 | `1431eb73c00e91f6b0ff893674daabc1b0090695c160d222454e0793bb898907` |
| ont/sample_ont_dupa_60_80/allele_2 | `f5b6acf1440ea35427c31bd7fcae251bbd8f6683a4cf1002b72edece25b6aceb` |
| ont/sample_ont_dupc_60_80/allele_1 | `1389fb8cb0b5001e2226f6d47845682978711cc78c4c37b6ff568b6d258fb3a7` |
| ont/sample_ont_dupc_60_80/allele_2 | `81e5aba3d33c1fe347f42badd76192b0f2415535a016d3a2cbcd7ff1a57f2298` |
| ont/sample_ont_normal_60_80/allele_1 | `372fb6f499ed4f4f07bb6274e5b2e8c1060bb6fcefec8311aecdef3b6d3661c8` |
| ont/sample_ont_normal_60_80/allele_2 | `118483bea7b612915fa140a4b8480e979552314d332ca979ab41873ce4d80a58` |
| heldout/sample_heldout_dupc_60_63/allele_1 | `a2a8c743bcb4e1129aea67d9e768ca5688a79ac09e06b6ebdda70806a473d2e6` |
| heldout/sample_heldout_dupc_60_63/allele_2 | `beb136233f14d032a140d90dde1c18d4aea90be461a1a746f3665e2cea2621fd` |
| heldout/sample_heldout_dupc_hap2_60_80/allele_1 | `d1dfcd7153fd26d3cc013b26110ca132b3265198c85a48a424c9b20c346c6e33` |
| heldout/sample_heldout_dupc_hap2_60_80/allele_2 | `9d8b1790c9abe7aa0a419b6fb20e413e5c995cf52c64f3222d6b7eecac106531` |
| heldout/sample_heldout_inscccc_terminal_60_80/allele_1 | `89f944847770b7307a9e66ba78d7e8c1736116a07b1431eefa6e073394367959` |
| heldout/sample_heldout_inscccc_terminal_60_80/allele_2 | `a2d73493e8485de3f5ca39df791efc5cfccd4a018f94b930edafdd2d94c21b95` |
| heldout/sample_heldout_normal_60_60/allele_1 | `7e5a963fabf490699e3a232445343d71cb3b60fba9261a81187e1adb3ab2a289` |
| heldout/sample_heldout_normal_60_60/allele_2 | `99c5aae4abe2a0ae3b4dc002cf55b5b203387b313525da4265c1dacaaacec551` |
| heldout/sample_heldout_normal_60_61/allele_1 | `5cc3581c8a5a4fdc62a51636ef4cf7348606930def3a2cf601bf84551948d71b` |
| heldout/sample_heldout_normal_60_61/allele_2 | `09b6104db34c02e7b967fd0fdc8a6886b6680b5bf868006ac959a584ad92b797` |
| heldout/sample_heldout_normal_60_62/allele_1 | `339e73acdcfc8541631e197a70a2e9451f2e77db4c65329ed3f1a88d5ae88123` |
| perturbations/sample_asymmetric_25_140__unique_names/allele_1 | `6d5e003f4697ab9e0335ca695990e9ce104549a45f760233056653583669c2b6` |
| perturbations/sample_asymmetric_25_140__unique_names/allele_2 | `7b84233007c700490393f0b37f7d3177863b202515bc98fc9f87a76f318b87d5` |
| perturbations/sample_dupc_100_120__unique_names/allele_1 | `3e6bd31e0d0433920e962c525ae5add366502dd2462c1629750716b3cdbffbdf` |
| perturbations/sample_dupc_100_120__unique_names/allele_2 | `c00f49bcd0ded877a4474be20aadaf7da039cdf6c865f1146d2950a6e5c9f1ea` |
| perturbations/sample_dupc_60_80__n20_seed1701/allele_1 | `ed341e050ba4afc66b9c0f7f28f43786e4c65e2cbe2d791f4f806cd446f7cf50` |
| perturbations/sample_dupc_60_80__n20_seed1702/allele_1 | `1389fb8cb0b5001e2226f6d47845682978711cc78c4c37b6ff568b6d258fb3a7` |
| perturbations/sample_dupc_60_80__n20_seed1703/allele_1 | `5e9acb5921b758d5cbe46be33d065334f6248da5d8ff1661b526826f707e4f84` |
| perturbations/sample_dupc_60_80__n20_seed1703/allele_2 | `99c5aae4abe2a0ae3b4dc002cf55b5b203387b313525da4265c1dacaaacec551` |
| perturbations/sample_dupc_60_80__n40_seed1701/allele_1 | `1389fb8cb0b5001e2226f6d47845682978711cc78c4c37b6ff568b6d258fb3a7` |
| perturbations/sample_dupc_60_80__n40_seed1701/allele_2 | `a583f234c8e51e58191ef7937a4a0d730535290f12bcf96acbe3d3fc0ff99955` |
| perturbations/sample_dupc_60_80__n40_seed1702/allele_1 | `1389fb8cb0b5001e2226f6d47845682978711cc78c4c37b6ff568b6d258fb3a7` |
| perturbations/sample_dupc_60_80__n40_seed1702/allele_2 | `05409fada27d913f87ca6b2d3cbb2b688160488469c98d8520a5070e161346f9` |
| perturbations/sample_dupc_60_80__n40_seed1703/allele_1 | `cfcbf00875c24f049507c3fed5c768b975e9d0659c31a92de97cafed56b20b5e` |
| perturbations/sample_dupc_60_80__n40_seed1703/allele_2 | `6f639b9b4d8947a4f9280c93cbf4f9e0d91e866a2ebdb92af5f190af55b11f17` |
| perturbations/sample_dupc_60_80__unique_names/allele_1 | `1389fb8cb0b5001e2226f6d47845682978711cc78c4c37b6ff568b6d258fb3a7` |
| perturbations/sample_dupc_60_80__unique_names/allele_2 | `57b1ed31030657ec26cdc3f4f91c0d8d346be01c1b6ee890c3951d4bb02d73e2` |
| perturbations/sample_homozygous_60_60__unique_names/allele_1 | `0bc9f181b1887836af8892ec35003f67531d5113d4b86710be9697717395d93e` |
| perturbations/sample_homozygous_60_60__unique_names/allele_2 | `87d23c2402d5d37ca4eac2131b61294e1cffcc0ec5b58a57252f77965746231e` |
| perturbations/sample_long_120_140__unique_names/allele_1 | `6207f5ebb01cff271dee8de53635cac4a6129a74256fd76a7896286fd2e8662e` |
| perturbations/sample_long_120_140__unique_names/allele_2 | `0d3d292beb08efed6a36480b2bac4f4d1f1c9c799567f8639b685e51116e1fe2` |
| perturbations/sample_normal_60_80__n20_seed1701/allele_1 | `937648394f3ae39435c1b71095b525fbd8ce0825ac1ba065bfc5cc6b8b1407b0` |
| perturbations/sample_normal_60_80__n20_seed1702/allele_1 | `5f74611229f79c47a35fbcfe0d143d5c625a174caa552f1fb014bc7bdef17628` |
| perturbations/sample_normal_60_80__n20_seed1702/allele_2 | `99c5aae4abe2a0ae3b4dc002cf55b5b203387b313525da4265c1dacaaacec551` |
| perturbations/sample_normal_60_80__n20_seed1703/allele_1 | `348c32c4f42d21fc92dff1829ceb556e6980b5a03496038428b0c6294c4063a2` |
| perturbations/sample_normal_60_80__n40_seed1701/allele_1 | `4b282e634e20fc4624c88374e39ff92860371b5eac5e9982d1353f1684b4bff3` |
| perturbations/sample_normal_60_80__n40_seed1701/allele_2 | `585134618eab09a882031c2d7aada653bba37222808b44415d09f52cc624a33e` |
| perturbations/sample_normal_60_80__n40_seed1702/allele_1 | `372fb6f499ed4f4f07bb6274e5b2e8c1060bb6fcefec8311aecdef3b6d3661c8` |
| perturbations/sample_normal_60_80__n40_seed1702/allele_2 | `99c5aae4abe2a0ae3b4dc002cf55b5b203387b313525da4265c1dacaaacec551` |
| perturbations/sample_normal_60_80__n40_seed1703/allele_1 | `bf964b212fdb9b2a2d356dcdc2917d34af9ad4e15985bb6bfa1dd5e602520f65` |
| perturbations/sample_normal_60_80__n40_seed1703/allele_2 | `d6c30b73ef63047b8e0b89b616035616587f9a703f2bd7aa6ef1899953126bfc` |
| perturbations/sample_normal_60_80__unique_names/allele_1 | `372fb6f499ed4f4f07bb6274e5b2e8c1060bb6fcefec8311aecdef3b6d3661c8` |
| perturbations/sample_normal_60_80__unique_names/allele_2 | `03f6a92decd602b6f254d1d71552ce6534dc262e7bd892112bc3b8012316722f` |
