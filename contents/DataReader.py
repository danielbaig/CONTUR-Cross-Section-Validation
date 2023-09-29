import os

import numpy as np


import re
from copy import deepcopy
import gzip



class DataReader:
    """
    A set of functions that read data from files.

    getHistogramLabel
    getSeed
    getParameters

    getBranchingRatios
        readBranchingRatios
        filterDecays
        addNewDecays
        doInvert
            substringInside
        checkDecays
        getOppositeCharges

    getCrossSections
        splitProductString
        crossSectionConvert
        uncertaintySplit

    getHistogramPoints

    getSubpool

    """

    def __init__(self, args):
        """
        Function run on class creation.

        Inputs:
            - args: Arguments when main.py was called.
        """

        
        self.pointIndex = args.point_index
        self.beam_energy = args.beam_energy
        self.run_area = args.run_area
        self.scan_area = args.scan_area

        # You will need to adjust this path for your needs.
        self.path = self.run_area + self.scan_area + f"/{self.beam_energy}TeV/{self.pointIndex}/"
        self.getParameters()

        if self.parameters['PionMass']>122 and self.parameters['PionMass']<222:
            tthad_selection = (('t','tbar'),('W+', 'b', 'bbar'), ('W-', 'b', 'bbar'), ('b','bbar'))
        else:   
            tthad_selection = (('t','tbar'),)

        
        possibleProductComb = {'TTHAD' : tthad_selection,
                            'MMJET' : (('Z0',),),
                            'L1L2METJET' : (('W+', 'W-'),),
        }

        
        
        try:
            self.getPool()
        except:
            raise Exception("Failed to get the pool information, ensure `contur runpoint_xxxx.yoda' "
                            + "and `contur-mkhtml --all --all-errbars' has been run inside the point's "
                            + "directory.")

        __,self.beam_energy,self.detection_mode = self.pool.split('_') # beam_energy in TeV
        self.desiredProducts = possibleProductComb[self.detection_mode]


        self.getHistogramLabel()
        self.readingSet = f"{self.histogramLabel.split('_')[0]}_{self.beam_energy}_{self.detection_mode}"

        self.getSeed()

        print(f'Looking in: {self.pool}/{self.histogramLabel}')
        

    def getPool(self):
        """
        Gets the most sensitive pool from the Summary.txt file.

        THIS SHOULD BE ALTERED TO NOT REQUIRE THE BEAM ENERGY AND OBTAIN THIS DIRECTLY FROM the .db FILE.
        """


        filePath = self.path + 'ANALYSIS/Summary.txt'

        i = 0
        line = True
        maxFileSize = int(1e+4)
        histLabels = []
        isReading = False

        possiblePools = []
        exclusions = []
        backgrounds = []
        

        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                if line[:5] == f'Pool:':
                    isReading = True
                    possiblePools.append(line[5:-1])
                    exclusions.append([])
                    backgrounds.append([])
                    i += 1
                    continue
                elif isReading and line == '\n':
                    # Ordered as DATABG, EXP, SMBG
                    # https://stackoverflow.com/questions/6618515/sorting-list-according-to-corresponding-values-from-a-parallel-list
                    exclusions[-1] = [x for _, x in sorted(zip(backgrounds[-1], exclusions[-1]))]
                    isReading = False
                    i += 1
                    continue
                
                # Read in the exclusion
                if isReading and line[:9]=='Exclusion':
                    thisLine = line.rstrip().split('=')
                    exclusions[-1].append(float(thisLine[1]))
                    backgrounds[-1].append(line[15:line.index('=')])
                # If there is no data, use an impossible entry.
                elif isReading and line[:2]=='No':
                    exclusions[-1].append(-1.)
                    backgrounds[-1].append(line[28:-1])


                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.') 


        # Using SM as bg
        exclusions = np.asarray(exclusions).T[2]

        maxIndex = np.argwhere(exclusions==np.max(exclusions))[0,0]
        self.pool = possiblePools[maxIndex]




    def getHistogramLabel(self):
        """
        Gets the histogram label from the Summary.txt file.
        """

        filePath = self.path + 'ANALYSIS/Summary.txt'

        i = 0
        line = True
        maxFileSize = int(1e+4)
        histLabels = []
        isReading = False
        

        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                if line[:-1] == f'Pool:{self.pool}':
                    isReading = True
                    i += 1
                    continue
                elif isReading and line == '\n':
                    isReading = False
                    i += 1
                    continue

                if isReading and line[0]=='/':
                    thisLine = line.rstrip().split('/')
                    histLabels.append(thisLine[1])


                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.') 

        if not np.all((hl==histLabels[0] for hl in histLabels)):
            raise Exception(f'Most sensitive histograms are not the same, this situation has not been implemented yet: {histLabels}')
        else:
            self.histogramLabel = histLabels[0]


    def getSeed(self):
        """
        Gets the seed number from the runpoint_xxxx.sh file.
        """

        filePath = self.path + f'runpoint_{self.pointIndex}.sh'

        i = 0
        line = True
        maxFileSize = int(1e+4)
    
    
        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                if line[:21] == 'Herwig run herwig.run':
                    start = line.index('=')
                    self.seed = line[start+1:start+4]
                    break

                i += 1


                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')


    def getParameters(self):
        """
        Get parameters of the model.
        """

        filePath = self.path + 'params.dat'
        i = 0

        line = True
        maxFileSize = int(1e+2)
        self.parameters = {}

        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                
                if len(line)!=0:
                    entries = line.split(' = ')
                    self.parameters[entries[0]] = float(entries[1][:-1])

                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')


    ################### Branching ratios ###################


    def getBranchingRatios(self):
        """
        Get branching ratios of decays that lead to any of the desired products.
    
        Outputs:
            - decays:list(2d): A nested list in format [[decayingParticle, product1, product2...],
                               partial width, branching ratio, isAllowed]
            - requiredDecayParticles:list: Particles that decay into the desired particles.
        """
    
        self.readBranchingRatios()
        

        # Add SM decays
        # https://www.physicsmasterclasses.org/exercises/keyhole/en/projects/number_of_families.html#:~:text=The%20official%20values%20of%20the,pairs%3A%20(3.360%2B%2F%2D%200.015)%25
        # OPAL Collaboration, G Abbiendi et al., ``W$ ^{+}$W$ ^{-}$ Production Cross section and W Branching Fractions in e$ ^{+}$e$ ^{-}$ collisions at 189 GeV,'' CERN-EP-2000-101, submitted to Phys. Lett. B.

        SM_decays = [[['Z0', 'mu+','mu-'], None, '0.03367', 'Yes'],
                    [['Z0', 'e+','e-'], None, '0.03366', 'Yes'],
                    [['Z0', 'tau+','tau-'], None, '0.03360', 'Yes'],
                    [['W+','e+','nu_e'], None, '0.1046', 'Yes'],
                    [['W+','mu+','nu_mu'], None, '0.1050', 'Yes'],
                    [['W+','tau+','nu_tau'], None, '0.1075', 'Yes'],]

        if self.detection_mode not in ('MMJET','TTHAD','L1L2METJET'):
            self.decays += SM_decays

 
        self.requiredDecayParticles = list(set([decay[0][0] for decay in self.decays]))
        

        self.filterDecays()
        self.getOppositeCharges()
        self.checkDecays()

        # Group particles that must decay from the same parent.
        if self.detection_mode == 'TTHAD' and (self.parameters['PionMass']>120 or self.parameters['PionMass']<220):
            for i,decay in enumerate(self.decays):
                if (len(decay[0][1:]) == 3 and ('W+' in decay[0][1:] or 'W-' in decay[0][1:]) 
                    and 'b' in decay[0][1:] and 'bbar' in decay[0][1:]):
                    self.decays[i][0] = [decay[0][0], 'Wbb']
                elif (len(decay[0][1:])==2 and 'b' in decay[0][1:] and 'bbar' in decay[0][1:]):
                    self.decays[i][0] = [decay[0][0], 'bbbar']

        elif self.detection_mode == 'L1L2METJET':
            for i,decay in enumerate(self.decays):
                if len(decay[0][1:])==2 and 'W+' in decay[0][1:] and 'W-' in decay[0][1:]:
                    decay[0] = [decay[0][0], 'W+W-']





    def readBranchingRatios(self):
        """
        Reads the branching ratios and decays from the .log file.
        """

        filePath = self.path + f'herwig-S{self.seed}-runpoint_{self.pointIndex}.log'

        i = 0
        isReading = False
        self.decays = []
        line = True
        maxFileSize = int(1e+4)
    
    
    
        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
    
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                if line[:9] == '# Parent:':
                    isReading = True
                    i += 2
                    __ = f.readline()
                    continue
                elif isReading and line[0] == '#':
                    isReading = False
                    i += 1
                    continue


                if isReading:
                    thisLine = line.rstrip().split(' ')
                    thisLine = list(filter(None,thisLine))
                    # https://stackoverflow.com/questions/4998629/split-string-with-multiple-delimiters-in-python
                    thisLine = [re.split(',|->',thisLine[0][:-1])] + thisLine[1:]
                    # https://stackoverflow.com/questions/740287/how-to-check-if-one-of-the-following-items-is-in-a-list
                    if thisLine[-1]=='Yes':
                        self.decays.append(thisLine)


                i += 1


                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')



    def filterDecays(self):
        """
        Only includes the decays which have non-zero cross sections and are necessary.
        """

        includedDecays = []
        print(self.desiredProducts)

        # Only include branching ratios that decay to the desired products or to another branching ratio
        for decay in self.decays:
            isNeeded  = (np.any([dp in decay[0][1:] for dp in self.requiredDecayParticles]) 
                            or np.any([np.any([dp in decay[0][1:] for dp in sub_desiredProducts])
                                                        for sub_desiredProducts in self.desiredProducts]))
  
            if isNeeded and float(decay[2]) > 0.:
                includedDecays.append(decay)

        self.decays = includedDecays


    def getOppositeCharges(self):
        """
        Obtains the opposite charged required particles and their decay paths.
        """

        new_rdps = []
        for rdp in self.requiredDecayParticles:
            if '+' in rdp:
                # Add oppositely charged particle.
                new_rdps.append(rdp.replace('+','-'))
        self.requiredDecayParticles += new_rdps

        # Get the length of all the particles to find
        charLen = np.fromiter(map(len, self.requiredDecayParticles), int)

        self.addNewDecays(new_rdps)


    def addNewDecays(self,new_rdps:list) -> list:
        """
        Adds the decays of opposite charges to the list of branching ratios.

        Inputs:
            - new_rdps:list: the new particle decays to add.
        """


        newDecays = []

        for new_rdp in new_rdps:
            # Search for original parity.
            lookFor = new_rdp.replace('-','+')
            for decay in self.decays:
                if decay[0][0] == lookFor:
                    newDecays.append(self.doInvert(decay))
                

        self.decays += newDecays

    
    def doInvert(self,decay:list) -> list:
        """
        Invert the decay to the oppositely charged decaying particle.

        Inputs:
            - decay:list: Decay information to invert.

        Outputs:
            - new_decay:list: Inverted decay.
        """

        new_decay = deepcopy(decay)
        neutralBSM = [rdp for rdp in self.requiredDecayParticles if '+' not in rdp and '-' not in rdp]

        for i in range(len(decay[0])):
            if '+' in decay[0][i]:
                new_decay[0][i] = decay[0][i].replace('+','-')
            elif '-' in decay[0][i]:
                new_decay[0][i] = decay[0][i].replace('+','-')
            elif self.substringInside('bar',decay[0][i]):
                new_decay[0][i] = decay[0][i][:-3]
            elif decay[0][i] in ('u','d','c','s','t','b','nu_e','nu_mu','nu_tau'):
                new_decay[0][i] = decay[0][i] + 'bar'
            elif decay[0][i] in neutralBSM + ['H']:
                pass
            else:
                raise Exception(f'Unrecognised product {decay[0][i]} for the decay: {decay}.')
                

        return new_decay

    @staticmethod
    def substringInside(substring:str, string:str):
        """
        Tests whether a substring is inside a longer string.
        
        Inputs:
            - substring:str: String that is being checked for.
            - string:str: String that contains the substring.

        Outputs:
            - isInside:bool: Whether the string is inside or not.
        """


        if len(string) < len(substring):
            return False


        for i in range(len(string) - len(substring)+1):
            if substring == string[i:i+len(substring)]:
                return True

        return False


    
    def checkDecays(self):
        """
        Ensure that the branching ratios add to one. If they do not, add an 'other' particle that
        makes up the difference.
        """

        threshold = 1e-6

        BR_sum = np.zeros(len(self.requiredDecayParticles))

        for i,rdp in enumerate(self.requiredDecayParticles):
            for decay in self.decays:
                if decay[0][0]==rdp:
                    BR_sum[i] += float(decay[2])

            # Threshold test
            if (abs(BR_sum[i] - 1.) > threshold and BR_sum[i] > 0.5
                or abs(BR_sum[i]) > threshold and BR_sum[i] < 0.5):
                ## A partial width for this may need to be implemented in the future.
                self.decays.append([[rdp, 'other'], None, 1 - BR_sum[i], 'Yes'])
               

    ################# Cross sections #####################
    
    def getCrossSections(self) -> list:
        """
        Get the cross section information from a file.

        Outputs:
            - cross_sections:list: The products of the collision and all non-zero
                                   cross sections in the format [products, cross sections].
        """

        self.readCrossSections()
        stringSplitter = lambda x: [self.splitProductString(x[0], self.requiredDecayParticles), x[1]]
        
        self.cross_sections = list(map(stringSplitter, self.cross_sections))

        # Remove uncertainties.
        self.crossSectionConvert()


    def readCrossSections(self):
        """
        Gets the cross section and collision product information from the .out file.
        """

        filePath = self.path + f'herwig-S{self.seed}-runpoint_{self.pointIndex}.out'

        i = 0
        isReading = False

        self.cross_sections = []
        maxFileSize = int(1e+4)
        line = True

        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                if line[:28] == 'Per matrix element breakdown':
                    isReading = True
                    i += 1
                    continue
                elif isReading and line[:5] == '=====':
                    isReading = False
                    i += 1
                    continue

                if isReading:
                    thisLine = line.rstrip().split(' ')
                    thisLine = list(filter(None,thisLine))
                    # Only include the products, the reactant partons are not required.
                    thisLine = [thisLine[0][thisLine[0].rindex('2')+1:], thisLine[-1]]
                    if thisLine[-1]!='0':
                        self.cross_sections.append(thisLine)


                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')



    @staticmethod
    def splitProductString(productString:str, BSM_particles:list):
        """
        Splits a product string into constituent decay products.

        Inputs:
            - productString:str: The string to break up.
            - BSM_particles:list: The BSM particles being added.

        """

        subsetSM = ['u','d','c','s','t','b','nu_e','nu_mu','nu_tau']
        SM_particles = (subsetSM 
                        + ['mu+','mu-','e+','e-','tau+','tau-'] 
                        + ['g','W+','W-','H','Z0','gamma']) 


        splits = [0]
        i = 0
        # Iterate forwards through the list.
        while i < len(productString):
            j = 1
            while j < len(productString) + 1 - i:
                if productString[i:i+j] in SM_particles+BSM_particles:
                    # Additionally check whether there is a bar afterwards too.
                    if productString[i+j:i+j+3] == 'bar':
                        barCheck = 3
                    # Check if g is actually the beginning of gamma.
                    elif productString[i+j:i+j+4] == 'amma':
                        barCheck = 4
                    else:
                        barCheck = 0

                    splits.append(i+j+barCheck)
                    i += barCheck + j - 1
                    break
                j += 1
            i += 1


        splitString = [productString[splits[i]:splits[i+1]] for i in range(len(splits)-1)]
        
        return splitString



    @staticmethod
    def uncertaintySplit(number:str):
        """
        Split a number in bracketed uncertainty format to a number and uncertainty separately.

        Inputs:
            - number:str: Number to split up.

        Outputs:
            - measurement:float: Actual reading.
            - uncertainty:float: Uncertainty of the reading.
        """

        order = number[number.index('e'):]
        measurement = number[:number.index('(')]
        uncertainty = number[number.index('(')+1:number.index(')')]

        lengthDifference = len(measurement) - len(uncertainty)

        # Format uncertainty
        if '.' in measurement:
            lengthDifference -= 1
        uncertainty = '0'*lengthDifference + uncertainty

        if '.' in measurement:
            uncertainty = uncertainty[0] + '.' + uncertainty[1:]

        # Combine with the order
        measurement = float(measurement + order)
        uncertainty = float(uncertainty + order)

        return measurement, uncertainty


    
    def crossSectionConvert(self):
        """
        Splits the cross sections and uncertainties.
        """

        converter = lambda x: [x[0], *self.uncertaintySplit(x[1])]

        self.cross_sections = list(map(converter, self.cross_sections))




    ####################### Histograms ##############################

    def getHistogramPoints(self):
        """
        Gets data, theory and SM points from a set of histograms.

        Outputs:
            - dataPoints:list: Data points for each histogram in format
                               (fileIndex, dataType, datapoint entry, datapoint property)
            - filenames:list: List of the histogram names.
            - yodaAreas:np.ndarray: Areas of histograms of simulated events.
        """

        dirPath = self.path + f"contur-plots/{self.readingSet}/{self.histogramLabel}"
        filePaths = []

        try:
            # https://www.geeksforgeeks.org/how-to-iterate-over-files-in-directory-using-python/
            for filename in os.scandir(dirPath):
                if filename.is_file() and filename.path[-4:]=='.dat':
                    filePaths.append(filename.path)
        except:
            raise Exception(f'Cannot find .dat files in {dirPath}.')


        self.filenames = np.asarray(np.fromiter(map(lambda x: x[-15:-4], filePaths), dtype='S16'), dtype=str)


        maxFileSize = int(1e+3)
        self.dataPoints = []
        self.additionalBinWidths = []
        # Read data from each file.
        for filePath in filePaths:
            self.readHistogram(filePath, maxFileSize)
        self.additionalBinWidths = np.asarray(self.additionalBinWidths)
        
        self.readSimulatedCrossSections()




    def readHistogram(self,filePath:str, maxFileSize:int):
        """
        Read histogram data from a particular table and store the data points.

        Inputs:
            - filePath:str: File to read from
            - maxFileSize:int: The maximum number of lines long a file can be.
        """

        line = True
        isReading = False
        self.dataPoints.append([[],[],[]])
        i = 0
        histCounter = 0
        additionalBinWidth = 1.


        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(filePath, 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                if line[:6] == '# xlow':
                    isReading = True
                    i += 1
                    continue
                elif isReading and line[:13] == '# END HISTO1D':
                    isReading = False
                    
                    i += 1
                    histCounter += 1
                    continue

                if isReading:
                    # Get data from the line.
                    thisLine = line.rstrip().split('\t')
                    thisLine = list(np.asarray(list(filter(None,thisLine)),dtype=float))
                    # Multiply the values and uncertainties by the additional bin width.
                    for j in range(2,5):
                        thisLine[j] *= additionalBinWidth
                    # What part to actually read.
                    self.dataPoints[-1][histCounter].append(thisLine)

                # Get additional bin from title.
                if line[:5] == 'Title' and '[' in line:
                    interval = line[line.index('[')+1:line.index(']')]
                    interval = interval.split(',')
                    additionalBinWidth = float(interval[1]) - float(interval[0])

                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')

        self.additionalBinWidths.append(additionalBinWidth)

   
    def readSimulatedCrossSections(self):
        """
        Read cross sections from the .yoda file. Requires getHistogramPoints is run first.
        """

        maxFileSize = int(1e+6)
        i = 0
        line = True
        yodaAreas = []

        try:
            thisFile = open(self.path + f'runpoint_{self.pointIndex}.yoda', 'r', encoding='UTF-8')
        except:
            thisFile = gzip.open(self.path + f'runpoint_{self.pointIndex}.yoda.gz', 'r')

        
        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with thisFile as f:
            while line and i < maxFileSize:
                line = f.readline().decode('utf-8')
                

                # Determines when to read from the document.
                if (line[:-13] == f'BEGIN YODA_HISTO1D_V2 /{self.histogramLabel}' 
                    and line[-12:-1] in self.filenames):
                    # Go to the relevant line.
                    for j in range(7):
                        i += 1
                        line = f.readline().decode('utf-8')
                    
                    if line[:7]=='# Area:':
                        area = float(line[8:-1])
                    else:
                        area = 0.


                    yodaAreas.append(area)
                    

                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')

        self.yodaAreas = np.asarray(yodaAreas)*self.additionalBinWidths


    
    def getSubpool(self) -> list:
        """
        Obtains the most sensitivie subpools from the Summary.txt file.

        Outputs:    
            - subpoolHist:list: The most sensitivie subpool.
        """

        self.subpoolHist = [[], [], []]


        i = 0
        counter = 0
        isReading = False

        maxFileSize = int(1e+4)
        line = True

        # https://stackoverflow.com/questions/3277503/how-to-read-a-file-line-by-line-into-a-list
        with open(self.path + 'ANALYSIS/Summary.txt', 'r', encoding='UTF-8') as f:
            while line and i < maxFileSize:
                line = f.readline()
                # Determines when to read from the document.
                # There is a \n at the end.
                if line[:-1] == f'Pool:{self.readingSet}':
                    isReading = True
                    i += 2
                    continue
                elif isReading and line == '\n':
                    isReading = False
                    i += 1
                    continue

                if isReading and line[0]=='/':
                    thisLine = re.split('/|\n|,', line)
                    # Only interested in the histogram names.
                    self.subpoolHist[counter] = thisLine[2::3]
                    counter += 1
                    
                i += 1

                if i >= maxFileSize:
                    raise Exception(f'Unsafe file size: number of lines exceeds {maxFileSize}.')




