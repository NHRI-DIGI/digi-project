mkdir /NovaSeq_128/Kenny/Digit_2024/PRS/PRS_score
cd PRS_score

Rscript /NovaSeq_128/Kenny/tool/bin/PRSice.R --dir . \
    --prsice /NovaSeq_128/Kenny/tool/bin/PRSice_linux \
    --base /NovaSeq_128/Kenny/Digit_2024/PRS/Basedata/Outputdata/base_completeqc.assoc.logistic \
    --target /NovaSeq_128/Kenny/Digit_2024/PRS/Targetdata/Outputdata/HBOC_p \
    --thread 1 \
    --stat OR \
    --binary-target T \
    --all-score
