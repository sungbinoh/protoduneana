import os

for i in range(11):
    ref_dedx = 1.5 + 0.1 * i
    ref_dedx_str = str(ref_dedx)
    print(ref_dedx_str)
    cp_fcl_template = "cp dEdx_calib_Birk_template.fcl dEdx_calib_Birk_" + ref_dedx_str + ".fcl"
    print(cp_fcl_template)
    os.system(cp_fcl_template)
    
    ## == Read normalization factors
    f_norm0 = open("../outputs/global_median_0_r5387_" + ref_dedx_str + ".txt")
    f_norm1 = open("../outputs/global_median_1_r5387_" + ref_dedx_str + ".txt")
    f_norm2 = open("../outputs/global_median_2_r5387_" + ref_dedx_str + ".txt")
    norm0 = ""
    norm1 = ""
    norm2 = ""
    for line in f_norm0:
        norm0 = line.split("\t")[1][:-1]
    for line in f_norm1:
        norm1 = line.split("\t")[1][:-1]
    for line in f_norm2:
        norm2 = line.split("\t")[1][:-1]

    print(norm0)
    print(norm1)
    print(norm2)

    ## == Write reference dE/dx (REFDEDX) and normalization factors (NORM0, NORM1, and NORM2)
    replace_REFDEDX = "find ./. -name dEdx_calib_Birk_" + ref_dedx_str + ".fcl -type f -exec sed -i s/REFDEDX/" + ref_dedx_str + "/g {} +"
    replace_norm0 = "find ./. -name dEdx_calib_Birk_" + ref_dedx_str + ".fcl -type f -exec sed -i s/NORM0/" + norm0 + "/g {} +"
    replace_norm1 = "find ./. -name dEdx_calib_Birk_" + ref_dedx_str + ".fcl -type f -exec sed -i s/NORM1/" + norm1 + "/g {} +"
    replace_norm2 = "find ./. -name dEdx_calib_Birk_" + ref_dedx_str + ".fcl -type f -exec sed -i s/NORM2/" + norm2 + "/g {} +"
    print(replace_REFDEDX)
    print(replace_norm0)
    print(replace_norm1)
    print(replace_norm2)
    os.system(replace_REFDEDX)
    os.system(replace_norm0)
    os.system(replace_norm1)
    os.system(replace_norm2)
