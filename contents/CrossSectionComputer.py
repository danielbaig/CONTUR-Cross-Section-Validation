
import numpy as np
from itertools import product,chain

from tqdm import tqdm



class CrossSectionComputer:
    """
    Computes cross sections from simulations and histograms.


    calculateCrossSection
        determineFinalParticles
            determineDecayProducts
            getCombinations
                makeDimensionsEven
                mergeProducts
        addToCrossSection
            getOccurances
    """
    
    def __init__(self, dataReader, isTesting:bool=False):
        """
        Initalises the class.

        Inputs:
            - cross_sections:list: All collision products with non-zero cross sections.
            - decays:list: Branching ratios for all relevant decays.
            - desiredProducts:tuple: Final products that are measured.
            - detection_mode:str: The method of particle detection.
            - isTesting:bool: Whether to run in testing mode.
            - parameters:dict: Parameters of the model.
        """

        self.cross_sections = dataReader.cross_sections
        self.decays = dataReader.decays
        self.desiredProducts = dataReader.desiredProducts
        self.detection_mode = dataReader.detection_mode
        self.isTesting = isTesting
        self.parameters = dataReader.parameters

        if self.detection_mode == 'TTHAD':
            if self.parameters['PionMass'] > 122 and self.parameters['PionMass']<222:
                self.desiredProducts = (('t','tbar'),('Wbb','Wbb'), ('Wbb','bbbar'),('bbbar','bbbar'))
            else:
                self.desiredProducts = (('t','tbar'),)
        elif self.detection_mode == 'L1L2METJET':
            self.desiredProducts = (('W+W-',),)

        # Determine the final branching ratio factors for each product type.
        _stopper = 1 if self.parameters['PionMass']>222 or self.parameters['PionMass']<122 else 4
        if self.detection_mode=='TTHAD':
            # https://www.hep.ucl.ac.uk/~jpc/all/ulthesis/node45.html            
            self.lastBranchingRatios = np.asarray((0.6832*0.6832, 0.6832*0.6832, 0.6832, 1.))[:_stopper]
        elif self.detection_mode=='MMJET':
            self.lastBranchingRatios = np.asarray((0.03367,))
        elif self.detection_mode=='L1L2METJET':
            self.lastBranchingRatios = np.asarray((2*0.1046*0.1050,))
        else:
            self.lastBranchingRatios = np.ones(len(self.desiredProducts))


        self.BR_U = 1e-4 ### Arbitrary magic number


    def calculateCrossSection(self):
        """
        Calculate the cross section for the given desired products.

        Outputs:
            - totalCrossSection:float: The calculated cross section for the given detection mode.
            - totalCrossSection_U:float: The uncertainty of the calculated cross section.
        """

        self.totalCrossSection = 0.
        self.totalCrossSection_U = 0.


        if self.isTesting:
            print('\nAll relevant decays:', self.decays)


        for cross_section in tqdm(self.cross_sections):
            finalParticles = self.determineFinalParticles(cross_section[0])
            self.addToTotalCrossSection(finalParticles, cross_section)

  
        self.totalCrossSection_U = np.sqrt(self.totalCrossSection_U)

        return self.totalCrossSection*1e+6, self.totalCrossSection_U*1e+6 # [fb]


    def determineFinalParticles(self,initialParticles:list) -> list:
        """
        Determine, using recursive methods, the final particles that are produced after decays.

        Note: The memory size of the list created is approximatly the same size as
              np.zeros(10,dtype=str).

        Inputs:
            - initialParticles:list: Particles produced from the partion collision.

        Outputs:
            - finalParticles:list: Final particles produced from the decays with their branching ratios.
        """

        try:
            # Finds all of the possible decay products.
            finalParticles = self.determineDecayProducts(initialParticles)
            
            if self.isTesting:
                print(f'\n======================= {initialParticles} ============================')
                print('\nStage 1:',finalParticles,'\n\n')
            # Determines all of the combinations of the decay products.
            finalParticles = self.getCombinations(finalParticles)


            # Tests if any relevant decays have occured.
            if [f for f in finalParticles if f!='+'] != initialParticles:
                # Combines the topmost layer.
                finalParticles, __ = self.makeDimensionsEven(finalParticles)


                # Combine all decaying particles' channnels.
                finalParticles = list(product(*finalParticles))
                finalParticles = list(map(lambda x: list(chain.from_iterable(x)), finalParticles))
            # In case no decays occured.
            else:
                finalParticles = [finalParticles]

            
            if self.isTesting:
                print('\nStage 2:',finalParticles,'\n\n')

        except:
            try:
                print(f'finalParticles: {finalParticles}')
            except:
                print('Unable to determine decay products.')


            raise Exception(f'Problem with {initialParticles}')

        return finalParticles



    def determineDecayProducts(self,inputParticles:list):
        """
        Determine the final particles and their branching ratios after all of the decays using recursion.

        Inputs:
            - inputParticles:list: A list of the particles to find the decay products of.

        Outputs:
            - finalProducts:list: The possible outcomes and their branching ratios.
        """

        finalProducts = []


        for inputParticle in inputParticles:
            decayFound = False
            # Search for available decays
            for decay in self.decays:
                if inputParticle == decay[0][0]:
                    if not decayFound:
                        finalProducts.append([])
                        decayFound = True

                    # Look at the next level of decay
                    output = self.determineDecayProducts(decay[0][1:])

                    finalProducts[-1].append(output)
                    # Add in the branching ratio
                    finalProducts[-1][-1].append('+' + str(decay[2]))


                
            if not decayFound:
                finalProducts.append(inputParticle)
      
        return finalProducts



    def getCombinations(self,nestlist:list) -> list:
        """
        Get all of the combinations of the final products for the set of decays.

        Inputs:
            - nestlist:list: The nested list to obtain the combinations from.

        Outputs:
            - nestlist:list: Combinations from the parent nestlist.
        """


        for i in range(len(nestlist)):
            if self.isTesting:
                print('item',nestlist[i])


            # Checks if the item is itself a string. s
            if (isinstance(nestlist[i],str)
                # Checks if the item only contains strings. los
                or np.all(np.fromiter(map(lambda x: isinstance(x,str), nestlist[i]),dtype=bool))
                # Check if the item is a list of lists of strings. lolos
                or (np.all(np.fromiter(map(lambda x: isinstance(x,list), nestlist[i]), dtype=bool)) 
                    and np.all([np.all([isinstance(item_sub,str) for item_sub in x]) for x in nestlist[i]]))):
                continue


            # Checks if the item contains any lists of lists. lol
            elif np.any([np.any([isinstance(x,list) for x in sub_item] for sub_item in nestlist[i])]):
                nestlist[i] = self.mergeProducts(nestlist[i])   

                if len(nestlist) >= 30:
                    raise Exception('List is of unsafe size.')

            else:
                raise Exception(f'An item in the list is not a list or a string:\n {nestlist[i]}')

        if self.isTesting:
            print('\nexiting',nestlist)

        
        return nestlist


    def mergeProducts(self, nestlist:list) -> list:
        """
        Merge a section of the nested list.

        Inputs:
            - nestlist:list: Part of the nested list.

        Outputs:
            - nestlist:list: Merged part of the nested list.
        """

        if self.isTesting:
            print('\nany lol')
        nestlist = self.getCombinations(nestlist)
        
        if self.isTesting:
            print('just bf',nestlist)

        
        # This ifelse chain uses the fact that combining and producting
        # the lists depends on the level at which they are being viewed
        # from in this particular method.
        
        if '+' in [nl[0] for nl in nestlist]:
            nestlist, dimsAdded = self.makeDimensionsEven(nestlist)

            nestlist = list(product(*nestlist))
            nestlist = list(map(lambda x: list(chain.from_iterable(x)), nestlist))

                       
            if self.isTesting:
                print('outcome',nestlist)

        else:
            nestlist, dimsAdded = self.makeDimensionsEven(nestlist)
            nestlist = list(chain(*nestlist))

            if self.isTesting:
                print('indicator outcome', nestlist)


        return nestlist


    @staticmethod
    def makeDimensionsEven(l:list):
        """
        Make the nested lists all of the same dimensionality.

        Inputs:
            - l:list: List to make dimensionally even.

        Outputs:
            - l:list: Processed list.
            - dimsAdded:int: Dimensions added.

        """

        shapes = np.empty(len(l),dtype=int)

        # Determines the depth of each sublist.
        for j,subItem in enumerate(l):
            temp = subItem
            k = 1
            
            # Check if it is already a string.
            if not isinstance(temp,list):
                shapes[j] = 0
                continue

            while k < 10:
                if not isinstance(temp[0],list):
                    shapes[j] = k
                    break
                else:
                    # Look into a deeper level.
                    temp = temp[0]

                k += 1
            else:
                raise Exception('List has been nested unreasonably deep (depth >= 10).')

     

        maxDim = np.max(shapes)
        # Determines if the depths are the same.
        if not np.all(shapes - maxDim == 0):
            for j,subItem in enumerate(l):
                for k in range(maxDim - shapes[j]):
                    # Adds dimension.
                    l[j] = [l[j]]
        else:
            return l,0

        return l, np.max(maxDim - shapes)



    def addToTotalCrossSection(self,finalParticles:list, cross_section:list):
        """
        Adds the cross section for this collision to the total cross section.

        Inputs:
            - finalParticle:list: Particles produced after decays.
            - cross_section:list: Cross section information in the form: 
                                  [[initialParticles], crossSection, crossSection_U]
        """


        totalBR = 0.
        totalBR_U = 0.

        for i,subFinalParticles in enumerate(finalParticles):
            # Test if all required products are in the finalParticles list
            if np.any([np.all([dp in subFinalParticles for dp in sub_desiredProducts]) 
                                                for sub_desiredProducts in self.desiredProducts]):
                try:
                    actualBR = np.product([float(item[1:]) 
                                                    for item in subFinalParticles if item[0]=='+'])
                    occurances = self.getOccurances(subFinalParticles)
                    dotProd = np.dot(occurances, self.lastBranchingRatios)
                    dot_U = np.linalg.norm(occurances*self.BR_U)


                    totalBR += actualBR * dotProd
                    actualBR_U = self.BR_U*dotProd
                    dotProd_U = actualBR*dot_U

                    totalBR_U += actualBR_U*actualBR_U + dotProd_U*dotProd_U

        

                except:
                    print('cross section', cross_section)
                    raise Exception('Problem with: ', subFinalParticles)


        totalBR_U = np.sqrt(totalBR_U)
        self.totalCrossSection += totalBR * float(cross_section[1])
        # Have assumed constant uncertainty to the branching ratios thus far.
        crossSection_U = totalBR * float(cross_section[2])
        BR_U = totalBR_U * float(cross_section[1])
        self.totalCrossSection_U += crossSection_U*crossSection_U + BR_U*BR_U



    def getOccurances(self, subFinalParticles:list) -> list:
        """
        Get the occurances of desired products in the list of final particles without double counting.

        Inputs:
            - subFinalParticle:list: Final decay products.

        Outputs:
            - occurances:list: Occurances of the desired products in the final particles.


        """
        occurances = []
            
        for i,sub_desiredProducts in enumerate(self.desiredProducts):
            occurances.append([0. for d in sub_desiredProducts])
            for desiredProduct in sub_desiredProducts:
                x = subFinalParticles.count(desiredProduct)
                y = sub_desiredProducts.count(desiredProduct)

                # Ensure no double counting.
                for index in np.argwhere(np.asarray(sub_desiredProducts)==desiredProduct):
                    occurances[i][index[0]] += x/y

        # Detected occurances is constrained by the particle that appears least, e.g., given
        # infinite tbar for TTHAD the occurannces is equivalent to the number of t.
        occurances = np.fromiter(map(min, occurances), dtype=float)


        if np.any([o!=round(o) for o in occurances]):
            raise Exception("occurances are not all integers", occurances)


        return occurances





    


def crudeApproximation(detection_mode:str, decays:list, cross_sections:list):
    """
    A crude approximation for the total cross sction to try to veryify the other method.

    Inputs:
        - detection_mode:str: Method of particle counting.
        - decays:list: Decays and their branching ratios.
        - cross_sections:list: Collision products and their cross sections.
    """

    DP0_ttbar = 0.
    rho0_UFO_ttbar = 0.
    rho_UFO_DPDP = 0.
    DP0_HZ0 = 0.
    DP_HW = 0.
    DP_Wbb = 0.
    DP0_bbbar = 0.
    H_bbbar = 0.
    DP0_H = 0.
    DP_H = 0.


    for decay in decays:
        if decay[0]==['DP0','t','tbar']:
            DP0_ttbar = float(decay[2])
        elif decay[0]==['rho0_UFO','t','tbar']:
            rho0_UFO_ttbar = float(decay[2])
        elif decay[0]==['rho+_UFO','DP0','DP+']:
            rho_UFO_DPDP = float(decay[2])
        elif decay[0]==['DP0', 'H', 'Z0']:
            DP0_HZ0 = float(decay[2])
        elif decay[0]==['DP+','H','W+']:
            DP_HW = float(decay[2])
        elif decay[0]==['DP+','Wbb']:
            DP_Wbb = float(decay[2])
        elif decay[0]==['DP0', 'bbar']:
            DP0_bbbar = float(decay[2])
        elif decay[0]==['H','bbbar']:
            H_bbbar = float(decay[2])
        if decay[0][0]=='DP0' and 'H' in decay[0][1:]:
            DP0_H += float(decay[2])
        elif decay[0][0]=='DP+' and 'H' in decay[0][1:]:
            DP_H += float(decay[2])

    print('DP0_ttbar, rho0_UFO_ttbar, rho_UFO_DPDP, DP0_HZ0, DP_HW, DP_Wbb\n', DP0_ttbar, rho0_UFO_ttbar, rho_UFO_DPDP, DP0_HZ0, DP_HW, DP_Wbb)

    tcs = 0.
    H_stuff = 0.00177*2+0.00147*3+0.0008*3+0.0004*3

    for cs in cross_sections:
        for part in cs[0]:
            if detection_mode=='TTHAD':
                if 'rho+_UFO' == part or 'rho-_UFO'==part:
                    tcs += float(cs[1])*(rho_UFO_DPDP*DP0_ttbar)
                elif 'rho0_UFO'==part:
                    tcs += float(cs[1])*rho0_UFO_ttbar
                elif 'DP0'==part:
                    tcs += float(cs[1])*(DP0_ttbar + DP0_H*H_bbbar)
                elif 'DP+'==part or 'DP-'==part:
                    tcs += float(cs[1])*(DP_Wbb+ DP_H*H_bbbar)
                elif 'H'==part:
                    tcs += float(cs[1])*H_bbbar

            elif detection_mode=='MMJET':
                if 'DP0'==part:
                    tcs += float(cs[1])*(DP0_HZ0*(1+H_stuff) + 0.00301366*H_stuff)
                elif 'DP+'==part or 'DP-'==part:
                    tcs += float(cs[1])*(DP_HW+0.00112)*H_stuff
                elif 'rho+_UFO'==part or 'rho-_UFO'==part:
                    tcs += float(cs[1])*((DP0_HZ0*(1+H_stuff) + 0.00301366*H_stuff) 
                                                        + (DP_HW+0.00112)*H_stuff)
                elif 'rho0_UFO'==part:
                    tcs += float(cs[1])*(DP_HW+0.00112)*H_stuff*2
                elif 'Z0'==part:
                    tcs += float(cs[1])

    if detection_mode=='TTHAD':
        tcs *= 0.6832
    elif detection_mode=='MMJET':
        tcs *= 0.03367

    print('Crude approximation:', tcs*1e+6)



    







