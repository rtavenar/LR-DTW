import os

__author__ = 'Romain Tavenard romain.tavenard[at]univ-rennes2.fr'

path = "../xp/ucr_all/"

fp_out = open("params.data", "w")
for ds_name in os.listdir(path):
    if os.path.isdir(os.path.join(path, ds_name)):
        fname_train = os.path.join(path, ds_name, "%s_TRAIN" % ds_name)
        fname_test = os.path.join(path, ds_name, "%s_TEST" % ds_name)
        if os.path.exists(fname_train) and os.path.exists(fname_test):
            for gamma in [0., .5, 1., 5., 10., 50.]:
                if gamma > 0:
                    for entropy_reg in [0, 1]:
                        fp_out.write("%s %s %f %d\n" % (fname_train, fname_test, gamma, entropy_reg))
                else:
                    fp_out.write("%s %s %f 0\n" % (fname_train, fname_test, gamma))
fp_out.close()
