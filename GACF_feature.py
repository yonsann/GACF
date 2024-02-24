import pandas as pd 
import numpy as np
import math
import warnings
import time

## Get GACF Feature

warnings.filterwarnings("ignore")

class Route2d:
    
    ## Horizontal Route
    
    def __init__(self, Data, N):
        self.data = Data # Standard_data
        self.alignment_chain = self.get_chain() # Alignment chain
        self.alignment_info = self.get_info() # Alignment information
        self.N = N # The number of terms of the series expansion
    
    def get_chain(self):
        # Get alignment chain
        # e.g. A--B--C--B--A  -> [A1, B2, C3, B4, A5]
        return [str(self.data['AlignmentType'].iloc[i]) + 
                str(self.data['RoadSectionOrder'].iloc[i]) for i in range(self.data.shape[0])]

    def get_info(self):
        # Get alignment information
        infor_dict = {}
        for i in range(self.data.shape[0]):
            road = str(self.data['AlignmentType'].iloc[i]) + str(self.data['RoadSectionOrder'].iloc[i])
            infor_dict[road] = {'Length': self.data['Length'].iloc[i], 
                                'CurvatureConstant': self.data['CurvatureConstant'].iloc[i], 
                                'DirectionMark': self.data['DirectionMark'].iloc[i]}
        return infor_dict

    def get_length(self, RoadID):
        # Get the length of road section "RoadID"
        return self.alignment_info[RoadID]['Length']

    def get_curvature(self, RoadID):
        # Get the curvature constant of road section "RoadID"
        return self.alignment_info[RoadID]['CurvatureConstant']

    def get_direction(self, RoadID):
        # Get the direction of road section "RoadID"
        return self.alignment_info[RoadID]['DirectionMark']

class Carfollowing2D(Route2d):
    
    # Temporary variable
    
    def get_variables(self, Infor1, Infor2):
        self.v1, self.T1, self.k1, self.e1, self.l1 = Infor1
        self.v2, self.T2, self.k2, self.e2, self.l2 = Infor2
        self.length1, self.length2 = self.alignment_info[self.T1]['Length'], self.alignment_info[self.T2]['Length']
        self.thou1, self.thou2 = self.alignment_info[self.T1]['CurvatureConstant'], self.alignment_info[self.T2]['CurvatureConstant']

        if (self.T1[0] == 'A') and (self.T2[0] == 'A'):
            variable1, variable2 = self.internalCF_AA()
        elif (self.T1[0] == 'B') and (self.T2[0] == 'B') and (self.alignment_info[self.T1]['DirectionMark'] == 'in'):
            variable1, variable2 = self.internalCF_BB1()
        elif (self.T1[0] == 'B') and (self.T2[0] == 'B') and (self.alignment_info[self.T1]['DirectionMark'] == 'out'):
            variable1, variable2 =  self.internalCF_BB2()
        elif (self.T1[0] == 'C') and (self.T2[0] == 'C'):
            variable1, variable2 =  self.internalCF_CC()
        elif (self.T1[0] == 'A') and (self.T2[0] == 'B'):
            variable1, variable2 =  self.adjacentCF_AB()
        elif (self.T1[0] == 'B') and (self.T2[0] == 'A'):
            variable1, variable2 =  self.adjacentCF_BA()
        elif (self.T1[0] == 'B') and (self.T2[0] == 'C'):
            variable1, variable2 =  self.adjacentCF_BC()
        elif (self.T1[0] == 'C') and (self.T2[0] == 'B'):
            variable1, variable2 =  self.adjacentCF_CB()
        else:
            variable1, variable2 =  [], []
        return variable1, variable2
    
    def internalCF_AA(self):
        # Scenario 1
        
        x1 = self.l1
        y1 = 0.0
        theta1 = 0.0
        alpha1 = 0.0

        x2 = self.l2
        y2 = 0.0
        theta2 = 0.0   
        alpha2 = 0.0
             
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def internalCF_BB1(self):
        # Scenario 2_1
        
        x1 = self.EulerSpiral_x(self.thou1, self.l1, self.length1, self.N)
        y1 = self.EulerSpiral_y(self.thou1, self.l1, self.length1, self.N)
        theta1 = self.EulerSpiral_angle(self.thou1, self.l1, self.length1)
        alpha1 = theta1

        x2 = self.EulerSpiral_x(self.thou2, self.l2, self.length2, self.N)
        y2 = self.EulerSpiral_y(self.thou2, self.l2, self.length2, self.N)
        theta2 = self.EulerSpiral_angle(self.thou2, self.l2, self.length2)
        alpha2 = theta2
        
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def internalCF_BB2(self):
        # Scenario 2_2

        x1 = self.EulerSpiral_x(self.thou1, self.e1, self.length1, self.N)
        y1 = self.EulerSpiral_y(self.thou1, self.e1, self.length1, self.N)
        theta1 = self.EulerSpiral_angle(self.thou1, self.e1, self.length1)
        alpha1 = math.pi - theta1

        x2 = self.EulerSpiral_x(self.thou2, self.e2, self.length2, self.N)
        y2 = self.EulerSpiral_y(self.thou2, self.e2, self.length2, self.N)
        theta2 = self.EulerSpiral_angle(self.thou2, self.e2, self.length2)
        alpha2 = math.pi - theta2
        
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def internalCF_CC(self):
        # Scenario 3

        x1 = 0.0
        y1 = 0.0
        theta1 = 0.0
        alpha1 = theta1

        x2 = math.sin((self.l2-self.l1)*self.thou1)/self.thou1 # or self.thou2
        y2 = (1-math.cos((self.l2-self.l1)*self.thou1))/self.thou1 # or self.thou2
        theta2 = (self.l2-self.l1)*self.thou1 # or self.thou2
        alpha2 = theta2
        
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def adjacentCF_AB(self):
        # Scenario 4

        x1 = -self.e1
        y1 = 0.0
        theta1 = 0.0
        alpha1 = theta1

        x2 = self.EulerSpiral_x(self.thou2, self.l2, self.length2, self.N)
        y2 = self.EulerSpiral_y(self.thou2, self.l2, self.length2, self.N)
        theta2 = self.EulerSpiral_angle(self.thou2, self.l2, self.length2)
        alpha2 = theta2
        
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def adjacentCF_BA(self):
        # Scenario 5

        x1 = self.EulerSpiral_x(self.thou1, self.e1, self.length1, self.N)
        y1 = self.EulerSpiral_y(self.thou1, self.e1, self.length1, self.N)
        theta1 = self.EulerSpiral_angle(self.thou1, self.e1, self.length1)
        alpha1 = math.pi-theta1

        x2 = -self.l2
        y2 = 0.0
        theta2 = 0
        alpha2 = math.pi - theta2
        
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def adjacentCF_BC(self):
        # Scenario 6
        
        xp = self.EulerSpiral_x(self.thou1, self.length1, self.length1, self.N)
        yp = self.EulerSpiral_y(self.thou1, self.length1, self.length1, self.N)
        thetap = self.EulerSpiral_angle(self.thou1, self.length1, self.length1)

        x1 = self.EulerSpiral_x(self.thou1, self.l1, self.length1, self.N)
        y1 = self.EulerSpiral_y(self.thou1, self.l1, self.length1, self.N)
        theta1 = self.EulerSpiral_angle(self.thou1, self.l1, self.length1)
        alpha1 = theta1

        theta2 = thetap + self.thou2*self.l2
        alpha2 = theta2  
        x2 = xp - math.sin(thetap)/self.thou2 + math.sin(theta2)/self.thou2
        y2 = yp + math.cos(thetap)/self.thou2 - math.cos(theta2)/self.thou2
        
        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]

    def adjacentCF_CB(self):
        # Scenario 7
        xp = self.EulerSpiral_x(self.thou2, self.length2, self.length2, self.N)
        yp = self.EulerSpiral_y(self.thou2, self.length2, self.length2, self.N)
        thetap = self.EulerSpiral_angle(self.thou2, self.length2, self.length2)

        theta1 = thetap + self.thou1*self.e1
        alpha1 = math.pi - theta1
        x1 = xp - math.sin(thetap)/self.thou1 + math.sin(theta1)/self.thou1
        y1 = yp + math.cos(thetap)/self.thou1 - math.cos(theta1)/self.thou1

        x2 = self.EulerSpiral_x(self.thou2, self.e2, self.length2, self.N)
        y2 = self.EulerSpiral_y(self.thou2, self.e2, self.length2, self.N)
        theta2 = self.EulerSpiral_angle(self.thou2, self.e2, self.length2)
        alpha2 = math.pi - theta2

        return [x1, y1, theta1, alpha1], [x2, y2, theta2, alpha2]
        

    def EulerSpiral_x(self, curvature_cons, s, length, N):
        # the x-coordinate on Euler Spiral
        x = 0.0
        for n in range(0, N):
            x += math.pow(curvature_cons,2*n) * math.pow(-1,n) * math.pow(s, 4*n+1) / ((4*n+1) * math.factorial(2*n) * math.pow(4,n) *math.pow(length, 2*n))
        return x

    def EulerSpiral_y(self, curvature_cons, s, length, N):
        # the y-coordinate on Euler Spiral
        y = 0.0
        for n in range(1, N+1):
            y += math.pow(curvature_cons,2*n-1) * math.pow(-1,n-1) * math.pow(s, 4*n-1) / ((4*n-1) * math.factorial(2*n-1) * math.pow(2,2*n-1) *math.pow(length, 2*n-1))
        return y
    
    def EulerSpiral_angle(self, curvature_cons, s, length):
        # the tangent angle on Euler Spiral
        theta = curvature_cons*math.pow(s,2)/(2*length)
        return theta
    
    def get_features(self, Infor1, Infor2):
        # Get GACF features

        # Variables
        variable1, variable2 = self.get_variables(Infor1, Infor2)
        x1, y1, theta1, alpha1 = variable1
        x2, y2, theta2, alpha2 = variable2
        v1, v2 = np.abs(self.v1), np.abs(self.v2)
        
        # GACF Feature 1: speed
        v = v1

        # GACF Feature 2: lateral relative space
        if theta1 == 0.0:
            dsx = np.abs(y2-y1)
        elif (theta1 > 0) and (theta1 < math.pi/2):
            dsx = np.abs(math.tan(theta1)*(x2-x1) - (y2-y1))/math.sqrt(math.pow(math.tan(theta1), 2)+1)
        elif theta1 == math.pi/2:
            dsx = np.abs(x2-x1)
        else:
            dsx = -1

        # GACF Feature 3: longitudinal relative space
        if theta1 == 0.0:
            dsy = np.abs(x2-x1)
        elif (theta1 > 0) and (theta1 < math.pi/2):
            dsy = np.abs(1/math.tan(theta1)*(x2-x1) + (y2-y1))/math.sqrt(math.pow(1/math.tan(theta1), 2)+1)
        elif theta1 == math.pi/2:
            dsy = np.abs(y2-y1)
        else:
            dsy = -1

        # GACF Feature 4: lateral relative speed
        dvx = v2*math.sin(alpha2-alpha1)

        # GACF Feature 5: longitudinal relative speed
        dvy = v2*math.cos(alpha2-alpha1) - v1
        
        return [v, dsx, dsy, dvx, dvy]