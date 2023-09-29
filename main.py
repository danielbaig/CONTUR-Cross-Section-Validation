import sys
from argparse import ArgumentParser

# from contents import *
from contents.HistogramComputer import HistogramComputer
from contents.DataReader import DataReader
from contents.CrossSectionComputer import CrossSectionComputer

from contents.CrossSectionComputer import crudeApproximation




def getArgs():
    """
    Defines the optional arguments.

    Outputs:
        - Arguments listed with --help.
    """

    parser = ArgumentParser()
    
    parser.add_argument('-i', '--point-index', help='Point to calculate the cross sections for.')
    parser.add_argument('-t','--isTesting',
                        help='Whether to run on predefined cross sections and collision products in test mode.',
                        default=False)
    parser.add_argument('-r', '--run-area', help='The run-area which dictates the model being analysed.')
    parser.add_argument('-b', '--beam-energy', help='COM energy of the collision.', default='13')
    parser.add_argument('-s', '--scan-area', help='Scan area name.', default='scan_all')

    args = parser.parse_args()

    if len(args.point_index) != 4:
        raise Exception('Point name must be of length four.')
    if args.run_area[-1] != '/':
        raise Exception(f"Run area path must end with '/' but ended with {args.run_area[-1]}")
    if '/' in args.scan_area[0]:
        raise Exception("Scan area must not contain any '/' (it must not start or end with '/' either).")
        

    return args   




def main(*args, **kwds):
    """
    Main function calls.
    """

    sys.setrecursionlimit(100)
    args = getArgs()


    dataReader = DataReader(args)


    # Obtain the cross section and branching ratio information from the files and convert to a usable format.
    try:
        print(f'Locating cross section and branching ratio infomation in: {dataReader.path}')
        
        dataReader.getBranchingRatios()
        dataReader.getCrossSections()
  

    except:
        raise Exception(f'Cannot find herwig-S(seed)-runpoint_{args.point_index}.out and/or herwig-S(seed)-runpoint_{args.point_index}.log and/or ANALYSIS/Summary.txt')

    # Display all non-zero cross sections that have relevant outcomes.
    if False:
        print(decays)
        for cs in cross_sections:
            print(cs)


    csc = CrossSectionComputer(dataReader, args.isTesting)

    # crudeApproximation(dataReader.detection_mode,decays,cross_sections)


    if args.isTesting:
        # For testing purposes.
        csc.cross_sections = [[['rho+_UFO','bbar'], '1e-6', '1e-7'],
                          [['rho0_UFO','g'], '1e-6', '1e-7'],]

        totalCrossSection,totalCrossSection_U = csc.calculateCrossSection()
        print(f'Trial cross section = ({totalCrossSection:0.5f}+/-{totalCrossSection_U:0.5f}) fb')

        return 0


    else:
        print('Calculating cross sections.')
        totalCrossSection, totalCrossSection_U = csc.calculateCrossSection()
        print(f'Total cross section: ({totalCrossSection:0.5f}+/-{totalCrossSection_U:0.5f}) fb')

    # Get cross sections from the histograms.
    dataReader.getHistogramPoints()
    dataReader.getSubpool()

    histComp = HistogramComputer(dataReader)

    histComp.computeHistogramCrossSections()
    histComp.plotDifferences(totalCrossSection, totalCrossSection_U)
    
    return 0



if __name__=='__main__':
    main()

