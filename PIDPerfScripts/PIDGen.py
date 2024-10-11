###############################################################################
# (c) Copyright 2021 CERN for the benefit of the LHCb Collaboration           #
#                                                                             #
# This software is distributed under the terms of the GNU General Public      #
# Licence version 3 (GPL Version 3), copied verbatim in the file "COPYING".   #
#                                                                             #
# In applying this licence, CERN does not waive the privileges and immunities #
# granted to it by virtue of its status as an Intergovernmental Organization  #
# or submit itself to any jurisdiction.                                       #
###############################################################################
from __future__ import print_function
from __future__ import division
from past.utils import old_div
import argparse, math, Meerkat
from PIDPerfScripts.PIDGenUtils import get_argparser, defaults, make_output_tree, \
    get_fill_objects
from ROOT import TFile, TH1F, gRandom, TTree, OneDimPhaseSpace, \
    CombinedPhaseSpace, BinnedDensity, Logger, addressof, std, \
    MyStruct
import PIDGenExpert.Run1.Config as ConfigRun1
import PIDGenExpert.Run2.Config as ConfigRun2
from math import sqrt, log


def main():
    '''Main function to run PIDGen.'''
    parser = get_argparser()
    parser.print_help()
    args = parser.parse_args()

    print(args)

    infilename = args.input
    intree = args.tree
    outfilename = args.output
    pidvar = args.pidvar
    ptvar = args.ptvar
    pvar = args.pvar
    etavar = args.etavar
    ntrvar = args.ntrvar
    minpid = args.lowerpid
    config = args.config
    dataset = args.dataset
    variant = args.var
    seed = args.seed
    ntrscale = args.ntrscale
    addcalibstat = args.calibstat

    if not infilename:
        print("Usage: PIDGen.py [options]")
        #  print "  For the usage example, look at pid_resample.sh file"
        print("  Available PID configs are: ")
        print("    For Run1 : ")
        for i in sorted(ConfigRun1.configs.keys()):
            print("      ", i)
        print("    For Run2 : ")
        for i in sorted(ConfigRun2.configs().keys()):
            print("      ", i)
        quit()
    run_pid_gen(
        infilename=infilename,
        intree=intree,
        outfilename=outfilename,
        pidvar=pidvar,
        ptvar=ptvar,
        pvar=pvar,
        etavar=etavar,
        ntrvar=ntrvar,
        minpid=minpid,
        config=config,
        dataset=dataset,
        variant=variant,
        seed=seed,
        ntrscale=ntrscale,
        addcalibstat=addcalibstat,
        noclone=args.noclone,
        outtree=args.outtree)


def run_pid_gen(infilename=None,
                intree=defaults['tree'],
                outfilename=defaults['output'],
                pidvar=defaults['pidvar'],
                ptvar=defaults['ptvar'],
                pvar=defaults['pvar'],
                etavar=defaults['etavar'],
                ntrvar=defaults['ntrvar'],
                minpid=defaults['lowerpid'],
                config=defaults['config'],
                dataset=defaults['dataset'],
                variant=defaults['var'],
                seed=defaults['seed'],
                ntrscale=defaults['ntrscale'],
                addcalibstat=defaults['calibstat'],
                noclone=defaults['noclone'],
                outtree=defaults['outtree']):
    '''Run PIDGen with the given config. intree can be a TTree or the name of
    the tree to be retrieved from the file named infilename. Similarly, outtree
    can either be the name of the tree to be saved to the file outfilename
    (default: the same name as intree), or a TTree to which the PIDGen branches
    will be added.'''
    if not infilename and not isinstance(intree, TTree):
        raise ValueError(
            "If intree isn't a TTree then you must specify infilename!")

    if variant == "default":
        variant = "distrib"  # to do: change this name in CreatePIDPdf

    year = None
    run = None
    try:
        year = dataset.split("_")[1]
    except:
        print(
            'Dataset format "%s" not recognized. Should be {MagUp,MagDown}_[Year]'
            % dataset)
        quit()
    if year in ["2011", "2012"]:
        run = 1
    elif year in ["2015", "2016", "2017", "2018"]:
        run = 2
    else:
        print('Data taking year "%s" not recognized' % year)
        quit()

    print(year, run, dataset)

    if run == 1:
        calibfilename = ConfigRun1.eosrootdir + "/" + config + "/" + "%s_%s.root" % (
            dataset, variant)
        transform_forward = ConfigRun1.configs[config]['transform_forward']
        transform_backward = ConfigRun1.configs[config]['transform_backward']
        configs = ConfigRun1.configs
    else:
        calibfilename = ConfigRun2.eosrootdir + "/" + config + "/" + "%s_%s.root" % (
            dataset, variant)
        configs = ConfigRun2.configs()
        if 'gamma' in list(configs[config].keys()):
            gamma = configs[config]['gamma']
            if gamma < 0:
                transform_forward = "(1.-(1.-x)**%f)" % abs(gamma)
                transform_backward = "(1.-(1.-x)**%f)" % (1. / abs(gamma))
            else:
                transform_forward = "((x)**%f)" % abs(gamma)
                transform_backward = "((x)**%f)" % (1. / abs(gamma))
        else:
            transform_forward = configs[config]['transform_forward']
            transform_backward = configs[config]['transform_backward']

    pidmin = 0.
    pidmax = 1.
    if 'limits' in configs[config]:
        pidmin = configs[config]['limits'][0]
        pidmax = configs[config]['limits'][1]
    if minpid == None:
        minpid = pidmin
    else:
        minpid = float(minpid)
        if minpid < pidmin: minpid = pidmin

    # Calculate the minimum PID variable to generate (after transformation)
    x = pidmin
    pidmin = eval(transform_forward)
    x = pidmax
    pidmax = eval(transform_forward)
    x = minpid
    minpid = eval(transform_forward)

    pid_phsp = OneDimPhaseSpace("PIDPhsp", pidmin, pidmax)
    mom_phsp = OneDimPhaseSpace("MomPhsp", 5.5, 9.5)
    eta_phsp = OneDimPhaseSpace("EtaPhsp", 1.5, 5.5)
    ntr_phsp = OneDimPhaseSpace("NtrPhsp", 3.0, 6.5)
    pidmom_phsp = CombinedPhaseSpace("PIDMomPhsp", pid_phsp, mom_phsp)
    pidmometa_phsp = CombinedPhaseSpace("PIDMomEtaPhsp", pidmom_phsp, eta_phsp)
    phsp = CombinedPhaseSpace("FullPhsp", pidmometa_phsp, ntr_phsp)

    infile = None
    if not isinstance(intree, TTree):
        infile = TFile.Open(infilename)
        tree = infile.Get(intree)
    else:
        tree = intree
    if not tree:
        print("Ntuple not found!")
        quit()

    kde = BinnedDensity("KDEPDF", phsp, calibfilename)
    nentries = tree.GetEntries()

    s = MyStruct()

    newtree, outfile = make_output_tree(tree, noclone, outfilename, outtree)

    branches = [newtree.Branch(pidvar, addressof(s, "newpid"), pidvar + "/D")]
    if addcalibstat:
        branches.append(
            newtree.Branch(pidvar + "_calibstat", addressof(s, "hint"),
                           pidvar + "_calibstat/D"))
    fillobjs = get_fill_objects(newtree, outtree, branches)

    if infile:
        infile.cd()

    h = TH1F("h", "h", 100, minpid, pidmax)
    if seed != None:
        gRandom.SetSeed(int(seed))
    else:
        gRandom.SetSeed(0)

    #x_code = compile("i.%s" % pidvar, '<string>', 'eval')
    if noclone:
        tree.SetBranchStatus('*', False)
        tree.SetBranchStatus(ptvar, True)
        tree.SetBranchStatus(ntrvar, True)
    pid_code = compile(transform_backward, '<string>', 'eval')
    log_pt_code = compile("log(i.%s)" % ptvar, '<string>', 'eval')
    if etavar == None:
        p_code = compile("i.%s" % pvar, '<string>', 'eval')
        pt_code = compile("i.%s" % ptvar, '<string>', 'eval')
        if noclone:
            tree.SetBranchStatus(pvar, True)
    else:
        eta_code = compile("i.%s" % etavar, '<string>', 'eval')
        if noclone:
            tree.SetBranchStatus(etavar, True)
    if ntrscale:
        log_ntracks_code = compile("log(float(i.%s)*%f)" % (ntrvar, ntrscale),
                                   '<string>', 'eval')
    else:
        log_ntracks_code = compile("log(float(i.%s))" % ntrvar, '<string>',
                                   'eval')

    Logger.setLogLevel(1)

    n = 0

    counter_nocalib = 0
    counter_lowcalib = 0

    for i in tree:
        point = std.vector('double')(4)
        point[0] = (pidmin + pidmax) / 2.
        point[1] = eval(log_pt_code)
        point[3] = eval(log_ntracks_code)
        if etavar == None:
            point[2] = -math.log(
                math.tan(math.asin(old_div(eval(pt_code), eval(p_code))) / 2.))
        else:
            point[2] = eval(eta_code)

        h.Reset()
        kde.slice(point, 0, h)
        hint = h.Integral()

        if hint > 0:
            x = h.GetRandom()
            if hint < 10:
                counter_lowcalib += 1
        else:
            x = minpid + (pidmax - minpid) * gRandom.Rndm()
            counter_nocalib += 1

        s.newpid = eval(pid_code)
        s.hint = hint

        for obj in fillobjs:
            obj.Fill()

        if (n % 1000 == 0):
            print("Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, PIDGen=%f, CalibStat=%f" % \
               (n, nentries, point[1], point[2], point[3], \
                s.newpid, s.hint))

        n += 1

    if noclone:
        tree.SetBranchStatus('*', True)

    # If adding branches to an existing output tree, make sure
    # the number of entries is set correctly for the tree.
    if isinstance(outtree, TTree):
        newtree.SetEntries(tree.GetEntries())

    if infile:
        infile.Close()
    if outfile:
        outfile.cd()
        newtree.Write()
        outfile.Close()

    print("PID Resampling finished.")
    print("  Total number of events processed:      ", n)
    print("  Events with no calibration:            ", counter_nocalib)
    print("  Eventw with low calib. stats (0<n<10): ", counter_lowcalib)


if __name__ == '__main__':
    main()
