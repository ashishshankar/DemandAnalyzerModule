# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 11:07:55 2016

@author: ashank3x
"""
from scipy import stats
import numpy
import numpy as np
import statsmodels.api as sm
from statsmodels.stats.stattools import durbin_watson
from scipy.stats import chisquare
import pandas as pd
import csv,sys,os,re,random
class demandAnalyzer(object):
    def __init__(self,path,datainsert):
        '''Constructor to initialize the bound object.'''
        self.path = path
        self.datainsert=datainsert
        self.data=[]
        self.data =self.readFromCSVFile()
        self.bp=None
        #self.x = numpy.array([i for i in range(1,31)])
        self.x =numpy.array([i for i in range(1,len(datainsert)+1)])
        self.y = numpy.array(self.data)
        self.slope, self.intercept, self.r_value, self.p_value, self.std_err = stats.linregress(self.x,self.y)
        
    def readFromCSVFile(self):
         if os.path.exists(self.path):
             os.remove(self.path)
         print(self.path) 
         out = open(self.path, 'w')
         datacsv =[]
         for i in self.datainsert:
              datacsv.append([i])
         for row in datacsv:
             for column in row:
                 out.write('%d' % column)
             out.write('\n')
         out.close()
         with open(self.path,'r') as csvfile:
             readData = csv.reader(csvfile)
             for completeData in readData:
                 for element in completeData:
                    self.data.append(int(element))
         return self.data
         
    def sac(self,x, k=1):
        """
        Sample autocorrelation (As used in statistics with normalization)
        http://en.wikipedia.org/wiki/Autocorrelation
        Parameters
        ----------
        x : 1d numpy array
            Signal
        k : int or list of ints
            Lags to calculate sample autocorrelation for
        Returns
        -------
        res : scalar or np array
            The sample autocorrelation. A scalar value if k is a scalar, and a
            numpy array if k is a interable.
        """
        try:
            res = []
            for ki in k:
                res.append(self.sac(x, ki))
            return np.array(res)
        except:
            pass
        mx = np.mean(x)
        if k==0:
            n = np.sum((x-mx)*(x-mx))
        else:
            n = np.sum((x[:-k]-mx)*(x[k:]-mx))
        d = len(x) * np.var(x)
        return n/d
    def boxpierce(self,x, lags, alpha=0.1):
        """
        The Box-Pierce test for determining if the data is independently distributed.
        Parameters
        ----------
        x : 1d numpy array
            Signal to test
        lags : int
            Number of lags being tested
           
        Returns
        -------
        Q : float
            Test statistic
        """
        n = len(x)
        Q = 0
        for k in range(1, lags+1):
            Q += (self.sac(x, k)**2)
        Q = n*Q
        return Q
    def boxpiercetestStatistic(self):
        self.bp=self.boxpierce(self.y, 7, alpha=0.1)
        print('Box-Pierce test statistic===={}'.format(self.bp))
        return self.bp
    def boxpierce_p_value(self):
        boxpiercePValue=(1-stats.chi2.cdf(self.bp ,7))
        print('Box-Pierce test p-value===={}'.format(boxpiercePValue))
        return boxpiercePValue
    def trendlineIntercept(self):
        print('The Trend line intercept:===={}'.format(self.intercept))
        return self.intercept
    def trendlineSlope(self):
        print('The Trend line intercept:===={}'.format(self.slope))
        return self.slope
    def slope_pValue(self):
        print('The Slope p-value===={}'.format(self.p_value))
        return self.p_value
    def zeroDemandfraction(self):
        zero_frac = 0
        for i in list(range(1,len(self.y))):
            if self.y[i]==0:
                zero_frac = zero_frac + 1
        zero_frac = int(zero_frac / len(self.y))
        print('The Zero demand fraction===={}'.format(zero_frac))
        return zero_frac
    def discreteIndicatorFlag(self):
        self.discreteIndicatorFlag=True
        for i in range(0,len(self.y)):
            if((round(self.y[i],0)) != (self.y[i])):
                self.discreteIndicatorFlag=False
        print('The Discrete Indicator Flag is===={}'.format(self.discreteIndicatorFlag))
        return self.discreteIndicatorFlag
        
    def piosondistributionMLE(self):
        tsMean = (sum(self.y)/len(self.y))
        mle = tsMean
        print('The Poisson distribution MLE is===={}'.format(mle))
        return mle
        
    def chi_square_test_degrees_Of_Freedom(self):
        dof =(len(set(self.y))-1)
        print('The Chi-square test degrees of freedom is===={}'.format(dof))
        return dof
        
    def chi_square_statistic_for_poison_distribution(self):
        demandValue =list(set(self.y))
        global p
        FOValue=[]
        for i in range(1,11):
        	if i != 9:
        		FOValue.append((list(self.y).count(i)))
        xfList=[]
        for i in list(zip(demandValue,FOValue)):
        	xfList.append(i[0]*i[1])
        p=0.00183980058
        PXValue= [0.0311,0.0793,0.1348,0.1719,0.1753,0.149,0.1086,0.0692,0.02]
        FOValuesum =sum(FOValue)
        FT = map(lambda x:x*FOValuesum,PXValue)
        FT=list(FT)
        chsqrt=chisquare(FOValue,FT,8)
        cropTheNumberFromChiSquare = (re.search('[\d.]+',str(chsqrt)).group()) #Used regular expression here to crop the value
        print('The Chi-square test statistic for Poisson distribution===={}'.format(cropTheNumberFromChiSquare))
        return cropTheNumberFromChiSquare
        
    def chi_square_statistic_for_poison_distribution_p_value(self):
        print('Chi-square test for Poisson distribution p-value:===={}'.format(p))
        return p
        
    def __del__(self):
        print("********************Instance for demandAnalyzer is deleted***********************")
        print("--Demand Analyzer module finished")
dataInsert=[]
for i in range(500):
    dataInsert.append(random.randrange(500))
#print(dataInsert)   Uncomment this line if u want the random number from 1 to 500
dataInsert=[5,6,2,10,2,1,7,8,7,3,1,6,4,2,10,7,5,5,4,8,1,10,3,6,4,1,4,7,6,8] 
filePath=(os.path.join(os.getcwd(),'ashish.csv'))      
obj = demandAnalyzer(filePath,dataInsert)  #User can give the file path name here

obj.trendlineIntercept()  #Calculating Trend line intercept:
obj.trendlineSlope()   #Calculating Trend line slope:
obj.slope_pValue()     #Calculating Box-Pierce test p-value:
obj.zeroDemandfraction()  #Calculating Zero demand fraction:
obj.boxpiercetestStatistic() #Calculating Box-Pierce test statistic:
obj.boxpierce_p_value()   # Calculating Box-Pierce p value
obj.discreteIndicatorFlag() #Calculating Discreate Indicator Flag
obj.piosondistributionMLE() #Calculating pioson distribution MLE Value
obj.chi_square_test_degrees_Of_Freedom()    #Calculating Chi-square test degrees of freedom
obj.chi_square_statistic_for_poison_distribution() # Calculating Chi-square test statistic for Poisson distribution
obj.chi_square_statistic_for_poison_distribution_p_value() #Calculating Chi-square test for Poisson distribution p-value

del(obj)   #Deleting the instance after completing of calculation of demand Analyzer