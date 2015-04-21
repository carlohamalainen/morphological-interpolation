import os.path

for i in range(0, 95 + 1):
    f = '/export/nif02/uqchamal/scratch_00init_morph/workflow_component_%d/morpho_%d/actual_final_merge_mincs/merged.mnc' % (i, i,)

    if not os.path.exists(f):
        print i
