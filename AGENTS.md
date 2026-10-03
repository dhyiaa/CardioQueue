# Codex project instructions

Read `PROJECT_HANDOFF_FOR_ARVIN.md` before making project-wide changes. Treat it as the navigation guide, not as a replacement for checking the cited source artifacts.

For manuscript work, use these sources in order:

1. `paper_MS_SI/Arvin Feedback.txt` for the requested editorial direction.
2. `paper_MS_SI/MS_BMC_MEDICAL_GENOMICS_FINAL.md` for the current main text.
3. `paper_MS_SI/SI_BMC_MEDICAL_GENOMICS_FINAL.md` for current supplementary text.
4. `paper_MS_SI/COMPLETION_AUDIT.md`, `paper_MS_SI/tables/`, and the cited result manifests for numerical verification.
5. `paper_MS_SI/Manuscript Writing skills.txt` for the local evidence and writing contract.

Never infer missing author, ethics, funding, competing-interest, controlled-access, repository, license, or contribution details. Preserve explicit verification markers until an author supplies the facts.

Never describe CardioQueue as an ACMG/AMP classifier, a calibrated clinical probability, a patient-priority score, prospective clinical validation, or a replacement for expert adjudication. Keep the distinction between the CardioQueue primary model, legacy source-held-out model, and temporal CardioQueue model.

Do not commit or expose participant-linked HiRO/CASPER WES/VERDICT records. Do not redistribute SHaRe-derived row-level data without documented permission. Do not redistribute FoldX binaries. Follow `.gitignore` and inspect staged files for secrets, private data, and files larger than GitHub's limits before every push.

