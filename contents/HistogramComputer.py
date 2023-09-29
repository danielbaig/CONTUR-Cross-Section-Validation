import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

class HistogramComputer:
    
    """
    A set of functions to compute cross sections from the histograms and display with the most
    sensivitive subpools.

    computeHistogramCrossSections
    plotDifferences
        plotHistogramCrossSections
        plotCrossSectionLabels
        createLegend

    """

    def __init__(self, dataReader):
        """
        Initalises the class.

        Inputs:
            - dataPoints:list: 4d list of the data points for all histograms for this point and 
                               detection mode for all three types of background.
            - subpoolHist:list: List of histograms in the most sensitive subpool.
            - filenames:np.ndarray: The names of each of the histograms.
            - histogramLabel:str: Name of the file holding the data.
        """


        self.dataPoints = dataReader.dataPoints
        self.subpoolHist = dataReader.subpoolHist
        self.filenames = dataReader.filenames
        self.histogramLabel = dataReader.histogramLabel
        self.yodaAreas = dataReader.yodaAreas
        
        self.numHist = len(self.filenames)



    def computeHistogramCrossSections(self):
        """
        Computes the cross section of each histogram.
        """

        print('Getting histogram areas.')

        # Calculations for the area of each bin.
        boxArea = lambda x: (x[1] - x[0])*x[2]
        boxArea_U = lambda x: (x[1] - x[0]) * np.asarray(x[3:])

        
        self.histCrossSections = np.empty((len(self.dataPoints), len(self.dataPoints[0])))
        self.histCrossSections_U = np.empty((*self.histCrossSections.shape,2))


        # Calculate the total area of the histogram for each type of data.
        for i,dp in enumerate(self.dataPoints):
            for j in range(3): # Only three types
                self.histCrossSections[i,j] = sum(map(boxArea,dp[j]))
                self.histCrossSections_U[i,j] = np.linalg.norm(np.asarray(list(map(boxArea_U,dp[j]))))




    def plotHistogramCrossSections(self, ax):
        """
        Plots the cross section differences of the histograms.

        Inputs:
            - ax: Instance of the plot.

        Outputs:
            - ax: Modified instance of the plot.
        """


        colours = ('orange','g')

        crossSections = np.zeros(2)
        crossSections_U = np.zeros((2,2,1))


        # Plot individual histograms.
        for i,filename in enumerate(self.filenames):
            d_BSM = self.histCrossSections[i,2] - self.histCrossSections[i,1]
            d_BSM_U = np.hypot(self.histCrossSections_U[i,2], self.histCrossSections_U[i,1])[:,None]

            colour = 'c'

            # Add to total if in subpool.
            for j in range(2):
                if filename in self.subpoolHist[2*j]:
                    crossSections[j] += d_BSM
                    crossSections_U[j] += d_BSM_U*d_BSM_U
                    colour = colours[j]
                    
            
            ax.scatter(d_BSM, i, c=colour, marker='.', zorder=2)
            ax.scatter(self.yodaAreas[i], i, c='y', marker='x', zorder=1)


        crossSections_U = np.sqrt(crossSections_U)

        # Most sensitive subpool cross sections
        for i in range(2):
            #ax.errorbar(crossSections[i], -2-i, xerr=crossSections_U[i],
            ax.scatter(crossSections[i], -2-i, 
                        c=colours[i], label=('databg','SMbg')[i],marker='.')

        return ax


    def plotCrossSectionLabels(self,ax, totalCrossSection:float, totalCrossSection_U:float):
        """
        Creates the labels for the plot.

        Inputs:
            - ax: Instance of the plot.
            - totalCrossSection:float: The calculated cross section addition due to BSM.
            - totalCrossSection_U:float: The uncertainty of the calculated cross section addition.

        Outputs:
            - ax: Modified instance of the plot.
        """

        
        # Additional lines
        ax.axvline(totalCrossSection, c='r', lw=0.5)
        ax.axvspan(totalCrossSection-totalCrossSection_U, totalCrossSection+totalCrossSection_U,
                    alpha=0.2, color='b')
        ax.axhline(-1,c='k', lw=0.5)
        ax.grid(axis='x')

        # Creates appropiate labels
        ax.set_title(self.histogramLabel)
        ax.set_xlabel('Cross section difference [fb]')
        # https://stackoverflow.com/questions/43673884/change-x-axis-ticks-to-custom-strings
        ax.set_yticks([-2.5] + list(range(self.numHist)))
        ax.set_yticklabels(['Subpool'] + list(self.filenames))

        # Set axes limits
        ax.set_ylim(-5, self.numHist+1)


        ax = self.createLegend(ax)
        
        return ax



    def createLegend(self,ax):
        """
        Creates a partly custom legend.

        Inputs:
            - ax: Instance of the plot.

        Outputs:
            - ax: Modified instance of the plot.
        """

        handles, labels = plt.gca().get_legend_handles_labels()

        # https://stackoverflow.com/questions/39500265/how-to-manually-create-a-legend
        new_label = 'calculated'
        new_patch = mpatches.Patch(color='r', label=new_label)
        handles.append(new_patch)
        labels.append(new_label)

        new_label = 'Yoda areas'
        new_line = mlines.Line2D([],[],color='y', marker='x', label=new_label, ls='')
        handles.append(new_line)
        labels.append(new_label)

        new_label = 'Not in subpool'
        new_line = mlines.Line2D([],[],color='c', marker='.', label=new_label,ls='')
        handles.append(new_line)
        labels.append(new_label)


        # https://stackoverflow.com/questions/13588920/stop-matplotlib-repeating-labels-in-legend
        new_label = dict(zip(labels,handles))


        ax.legend(new_label.values(), new_label.keys(), loc='upper right')

        return ax




    def plotDifferences(self,totalCrossSection:float, totalCrossSection_U:float):
        """
        Displays the cross section measurements and compares the differences with the calculated ones.

        Inputs:
            - totalCrossSection:float: The calculated cross section addition due to BSM.
            - totalCrossSection_U:float: The uncertainty of the calculated cross section addition.
        """

       
        fig = plt.figure(figsize=(6,max(6,0.2*self.numHist)),dpi=200)
        ax = fig.add_subplot()

       
        ax = self.plotHistogramCrossSections(ax)
        ax = self.plotCrossSectionLabels(ax, totalCrossSection, totalCrossSection_U)
            
        
        histPath = '/home/dbaig/Documents/cross_section_calcs/histogram_comparisons.png'
        # https://stackoverflow.com/questions/9622163/save-plot-to-image-file-instead-of-displaying-it
        plt.savefig(histPath, bbox_inches='tight')

        print(f'Figure saved to: {histPath}')

