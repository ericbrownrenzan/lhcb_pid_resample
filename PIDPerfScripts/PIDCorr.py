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


from past.utils import old_div
import sys, math, Meerkat
from PIDPerfScripts.PIDGenUtils import get_argparser, defaults, make_output_tree, \
    get_fill_objects
from ROOT import TFile, TH1F, OneDimPhaseSpace, CombinedPhaseSpace, \
    BinnedDensity, Logger, TTree, MyStruct, addressof, std
import PIDGenExpert.Run1.Config as ConfigRun1
import PIDGenExpert.Run2.Config as ConfigRun2
import PIDGenExpert.Run1.ConfigMC as ConfigMCSim08
import PIDGenExpert.Run1.ConfigMCSim09 as ConfigMCSim09
import PIDGenExpert.Run2.ConfigMC as ConfigMCRun2
from math import sqrt, log


def main():
    '''Main function to run PIDCorr.'''
    parser = get_argparser(False)
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
    oldpidvar = args.simpidvar
    conf = args.config
    dataset = args.dataset
    variant = args.var
    simversion = args.simversion
    ntrscale = args.ntrscale
    addcalibstat = args.calibstat

    if not infilename:
        print("Usage: PIDCorr.py [options]")
        #  print "  For the usage example, look at pid_transform.sh file"
        print("  Available PID configs for Run1/sim08 are: ")
        for i in sorted(ConfigRun1.configs.keys()):
            if i in list(ConfigMCSim08.configs.keys()):
                print("    ", i)
        print("  Available PID configs for Run1/sim09 are: ")
        for i in sorted(ConfigRun1.configs.keys()):
            if i in list(ConfigMCSim09.configs.keys()):
                print("    ", i)
        print("  Available PID configs for Run2/sim09 are: ")
        for i in sorted(ConfigRun2.configs().keys()):
            if i in list(ConfigMCRun2.configs.keys()):
                print("    ", i)

        # Exit politely
        sys.exit(0)
    run_pid_corr(
        infilename=infilename,
        intree=intree,
        outfilename=outfilename,
        pidvar=pidvar,
        ptvar=ptvar,
        pvar=pvar,
        etavar=etavar,
        ntrvar=ntrvar,
        minpid=minpid,
        oldpidvar=oldpidvar,
        config=conf,
        dataset=dataset,
        variant=variant,
        simversion=simversion,
        ntrscale=ntrscale,
        addcalibstat=addcalibstat,
        noclone=args.noclone,
        outtree=args.outtree)


def run_pid_corr(infilename=None,
                 intree=defaults['tree'],
                 outfilename=defaults['output'],
                 pidvar=defaults['pidvar'],
                 ptvar=defaults['ptvar'],
                 pvar=defaults['pvar'],
                 etavar=defaults['etavar'],
                 ntrvar=defaults['ntrvar'],
                 minpid=defaults['lowerpid'],
                 oldpidvar=defaults['simpidvar'],
                 config=defaults['config'],
                 dataset=defaults['dataset'],
                 variant=defaults['var'],
                 simversion=defaults['simversion'],
                 ntrscale=defaults['ntrscale'],
                 addcalibstat=defaults['calibstat'],
                 noclone=defaults['noclone'],
                 outtree=defaults['outtree']):
    '''Run PIDCorr with the given config. intree can be a TTree or the name of
    the tree to be retrieved from the file named infilename. Similarly, outtree
    can either be the name of the tree to be saved to the file outfilename
    (default: the same name as intree), or a TTree to which the PIDGen branches
    will be added.'''
    if not infilename and not isinstance(intree, TTree):
        raise ValueError(
            "If intree isn't a TTree then you must specify infilename!")

    if simversion == "sim08":
        ConfigMC = ConfigMCSim08
        Config = ConfigRun1
    elif simversion == "sim09":
        ConfigMC = ConfigMCSim09
        Config = ConfigRun1
    elif simversion == "run2":
        ConfigMC = ConfigMCRun2
        Config = ConfigRun2
    else:
        print("Simulation version %s unknown" % simversion)
        sys.exit(1)

    if variant == "default":
        variant = "distrib"  # to do: change this name in CreatePIDPdf

    datapdf = Config.eosrootdir + "/" + config + "/" + dataset + "_" + variant + ".root"
    simpdf = ConfigMC.eosrootdir + "/" + config + "/" + dataset + "_" + variant + ".root"

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
            elif gamma == 1.:
                transform_forward = "x"
                transform_backward = "x"
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

    if isinstance(intree, TTree):
        tree = intree
        infile = None
    else:
        infile = TFile.Open(infilename)
        tree = infile.Get(intree)
    if not tree:
        print("Ntuple not found!")
        sys.exit(1)

    nentries = tree.GetEntries()

    datakde = BinnedDensity("KDEPDF", phsp, datapdf)
    simkde = BinnedDensity("KDEPDF", phsp, simpdf)

    s = MyStruct()

    newtree, outfile = make_output_tree(tree, noclone, outfilename, outtree)

    branches = [newtree.Branch(pidvar, addressof(s, "newpid"), pidvar + "/D")]
    if addcalibstat:
        branches.append(
            newtree.Branch(pidvar + "_calibstat", addressof(s, "hint"),
                           pidvar + "_calibstat/D"))
        branches.append(
            newtree.Branch(pidvar + "_mcstat", addressof(s, "hintmc"),
                           pidvar + "_mcstat/D"))
    fillobjs = get_fill_objects(newtree, outtree, branches)

    if infile:
        infile.cd()

    hdata = TH1F("hdata", "h", 100, minpid, pidmax)
    hsim = TH1F("hsim", "h", 100, minpid, pidmax)

    print(transform_backward)

    if noclone:
        tree.SetBranchStatus('*', False)
        tree.SetBranchStatus(oldpidvar, True)
        tree.SetBranchStatus(ptvar, True)
        tree.SetBranchStatus(ntrvar, True)
    var_code = compile("i.%s" % oldpidvar, '<string>', 'eval')
    pid_code = compile(transform_backward, '<string>', 'eval')
    oldpid_code = compile(transform_forward, '<string>', 'eval')
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
    counter_nomc = 0
    counter_lowcalib = 0
    counter_lowmc = 0

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

    #  print point[0], point[1], point[2], point[3]

        hdata.Reset()
        hsim.Reset()
        datakde.slice(point, 0, hdata)
        simkde.slice(point, 0, hsim)

        s.hint = hdata.Integral()
        s.hintmc = hsim.Integral()
        if s.hint == 0:
            counter_nocalib += 1
        elif s.hint < 10:
            counter_lowcalib += 1
        if s.hintmc == 0:
            counter_nomc += 1
        elif s.hintmc < 10:
            counter_lowmc += 1

        x = eval(var_code)
        oldpid = x
        if transform_forward == "x" or x >= 0:
            oldpid = eval(oldpid_code)
            if oldpid < pidmin or oldpid > pidmax:
                x = oldpid
            else:
                x = datakde.transform(hsim, hdata, oldpid)
            s.newpid = eval(pid_code)
        else:  # The case for ProbNN<0, just leave as it is
            s.newpid = x

        for obj in fillobjs:
            obj.Fill()

        if (n % 1000 == 0):
            print("Event %d/%d : Pt=%f, Eta=%f, Ntr=%f, OldPID=%f, PIDCorr=%f, X=%f, CalibStat=%f, MCStat=%f" % \
               (n, nentries, point[1], point[2], point[3], \
                oldpid, s.newpid, x, s.hint, s.hintmc))

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
    print("  Events with no MC:                     ", counter_nomc)
    print("  Eventw with low MC stats (0<n<10):     ", counter_lowmc)


if __name__ == '__main__':
    main()
