'''
Created on Jun 1, 2016

@author: ficar
'''

#mRow,mCol,mVal
#print str(list(qq))
#print matrix(qq)
#print matrix(numpy.diag(numpy.diag(qq)))
#print str(qq.echelon_form())
#matrix(numpy.diag(numpy.diag(qq))

'''Coded by Sal Ficara Version 2.0.0'''

import numpy
import math
import scipy
Pi = 3.141592

class SteepestDescents():


    def __init__(self, aList):
        self.itZero = aList
        self.itNumber = 0
        self.alpha3 = 1
        self.tolerance = .000001
        self.g3Val = None

    def __repr__(self):
        return "Instance of the Steepest Descents Class"

    def theJacobian(self, x1, x2, x3):
        if x3 == None:
            jacTwoEq = numpy.matrix(2,2,[6*x1, -2*x2, 3*x2^2-3*x1^2, 6*x1*x2])
            return jacTwoEq
        else:
            jacEqThreeEq = numpy.matrix(3,3,[3, x3*math.sin(x2*x3),0,2*x1,-162*x2-16.2, math.cos(x3),-x2*math.exp(-x1*x2),
                                       -x1*math.exp(-x1*x2),20])
            return jacEqThreeEq



    #For systems of TWO EQUATIONS enter equation 1 for the variable eq1 AND
    #enter equation 2 for the variable eq2
    def f1(self,x1, x2):

        eq1 = 3*x1^2 - x2^2
        eq2 = 3*x1*x2^2 - x1^3 -1

        g1 = eq1^2 + eq2^2
        return eq1, eq2, g1

    #for systems of THREE EQUATIONS enter equation 1 for the variable eq1, enter equation 2
    #for the variable eq2 and enter equation 3 or the variable eq3
    def f1Version2(self,x1,x2,x3):
        eq1 = 3*x1-math.cos(x2*x3) - 1/2
        eq2 = x1^2 - 81*(x2 + 0.1)^2 + math.sin(x3) + 1.06
        eq3 = math.e^(-x1*x2) + 20*x3 + (10*Pi - 3)/3
        g1 = eq1^2 + eq2^2 + eq3^2
        return eq1, eq2, eq3, g1



    def zFind(self,jac, xO1, xO2, xO3):
        firstMatrix = jac
        if xO3 == None:
            values = self.f1(xO1,xO2)
            Fx = numpy.matrix(2, 1,[values[0], values[1]])

            squareRN = 2*firstMatrix * Fx
            v1, v2 = squareRN[0,0], squareRN[1,0]
            zO = v1^2 + v2^2
            print("zO: " + str(zO))
            zOFinal = float(1/math.sqrt(zO))
            print ("\nz1: " + str(zOFinal*v1))
            print ("z2: " + str(zOFinal*v2))
            return round(zOFinal*v1,6), round(zOFinal*v2,6), None, values[2], xO1, xO2, None
        else:
            values = self.f1Version2(xO1,xO2,xO3)
            Fx = numpy.matrix(3, 1,[values[0], values[1], values[2]])

            squareRN = 2*firstMatrix * Fx
            v1, v2, v3 = squareRN[0,0], squareRN[1,0], squareRN[2,0]
            #gradVector = numpy.vector([v1,v2,v3])
            print ("\nGradient Vector: " + str(v1 , v2, v3))
            zO = v1^2 + v2^2 + v3^2
            zOFinal = float(1/math.sqrt(zO))
            print ("\nMake Gradient A Unit Vector" + str(zOFinal*v1,zOFinal*v2,zOFinal*v3))
            print ("\nz1: " + str(round(zOFinal*v1,6)))
            print ("z2: " + str(round(zOFinal*v2,6)))
            print ("z3: " + str(round(zOFinal*v3,6)))
            return round(zOFinal*v1,6), round(zOFinal*v2,6), round(zOFinal*v3,6), values[3], xO1, xO2, xO3

    def g3(self, alpha, tolerance, zVal, g1, initValue):
        timesDivided = 0
        if zVal.nrows() == 2:
            newX1 = initValue[0] - alpha*zVal[0,0]
            newX2 = initValue[1] - alpha*zVal[1,0]
            g3Val = self.f1(newX1, newX2)
            self.g3Val = g3Val[2]
            if self.g3Val < g1:
                print ("\n" + str(g3Val[2]) + " is less than " + str(g1)+ "\n")
            elif self.g3Val >= g1:
                print ("\n***ERROR***\n")
                print (str(g3Val[2]) + " greater than " + str(g1) + ". "+ "We need g3 value to be less than g1" "")
                print ("Taking corrective measures... ")
                while self.g3Val > g1:
                    self.alpha3 = self.alpha3/2

                    newX1 = initValue[0] - self.alpha3*zVal[0,0]
                    newX2 = initValue[1] - self.alpha3*zVal[1,0]
                    g3Val = self.f1(newX1, newX2)
                    self.g3Val = g3Val[2]

                    timesDivided += 1

                print ("ERROR CORRECTION: " + str(g3Val[3]) + " is less than " + str(g1))
                print ("TIMES DIVIDED: " + str(timesDivided))

            newX1 = initValue[0] - (self.alpha3/2)*zVal[0,0]
            newX2 = initValue[1] - (self.alpha3/2)*zVal[1,0]

            g2Val = self.f1(newX1, newX2)

            print ("\ng1: " + str(g1) + " W/ Alpha 0 = 0")
            print ("g2: " + str(g2Val[2]) + " W/ Alpha_2 = " + str(float(self.alpha3/2)))
            print ("g3: " + str(self.g3Val) + " W/ Alpha_3 = " + str(float(self.alpha3)))

            return g1, g2Val[2], self.g3Val, self.alpha3, self.alpha3/2
        else:

            newX1 = initValue[0] - alpha*zVal[0,0]
            newX2 = initValue[1] - alpha*zVal[1,0]
            newX3 = initValue[2] - alpha*zVal[2,0]
            g3Val = self.f1Version2(newX1, newX2, newX3)
            self.g3Val = g3Val[3]
            if self.g3Val < g1:
                print ("" + str(self.g3Val) + " is less than " + str(g1)+ "\n")
            elif self.g3Val >= g1:
                print ("\n***ERROR***\n")
                print (str(self.g3Val) + " greater than " + str(g1) + ". "+ "We need g3 value to be less than g1" "")
                print ("Taking corrective measures... ")
                while self.g3Val > g1:
                    self.alpha3 = self.alpha3/2

                    newX1 = initValue[0] - self.alpha3*zVal[0,0]
                    newX2 = initValue[1] - self.alpha3*zVal[1,0]
                    newX3 = initValue[2] - self.alpha3*zVal[2,0]
                    g3Val = self.f1Version2(newX1, newX2, newX3)
                    self.g3Val = g3Val[3]

                    timesDivided += 1

                print ("ERROR CORRECTION: " + str(g3Val[3]) + " is less than " + str(g1))
                print ("TIMES DIVIDED: " + str(timesDivided))



            newX1 = initValue[0] - (self.alpha3/2)*zVal[0,0]
            newX2 = initValue[1] - (self.alpha3/2)*zVal[1,0]
            newX3 = initValue[2] - (self.alpha3/2)*zVal[2,0]

            g2Val = self.f1Version2(newX1, newX2, newX3)

            print ("\ng1: " + str(g1) + " W/ Alpha 0 = 0")
            print ("g2: " + str(g2Val[3]) + " W/ Alpha_2 = " + str(float(self.alpha3/2)))
            print ("g3: " + str(g3Val[3]) + " W/ Alpha_3 = " + str(float(self.alpha3)))

            return g1, g2Val[3], self.g3Val, self.alpha3, self.alpha3/2



    def findHVal(self,g1, g2, g3, alpha3, alpha2):
        h1 = (g2-g1)/(alpha2 - 0)
        h2 = (g3-g2)/(alpha3 - alpha2)
        h3 = (h2 - h1)/(alpha3 - 0)
        print ("\nh1 : " + str(h1))
        print ("h2 : " + str(h2))
        print ("h3 : " + str(h3))
        return h1, h2, h3, alpha2, g1

    def createNewtonIntPoly(self,h1, h2, h3, alpha2, g1):
        p1 = numpy.poly1d([g1])
        p2 = numpy.poly1d([h1, 0])
        p3 = numpy.poly1d([h3,0])
        p4 = numpy.poly1d([1 , -alpha2])
        p5 = numpy.polymul(p3,p4)
        finalTerm = p1 + p2 + p5
        finalTermAltWay =lambda x: g1 + h1*x + h3*x*(x-alpha2)
        
        #print ("The MIN" + str(find_local_minimum(finalTermAltWay, 0, self.alpha3)))
        #need to make a function for finding minimum
        b = scipy.optimize.local(finalTermAltWay, 0, self.alpha3)

        print ("\nNewtonian Int. Root: " + str(float(numpy.roots(numpy.polyder(finalTerm)))))
        #return b[1]
        return numpy.roots(numpy.polyder(finalTerm))

    def findFirstIte(self,root, xO1, xO2, xO3, z1, z2, z3):
        if z3==None:
            #print root, xO1, xO2, z1, z2
            x1 = xO1 - root*z1
            x2 = xO2 - root*z2
            print ("\n")
            print ("x1: ",float(x1))
            print ("x2: ",float(x2))
            return float(x1),float(x2), None

        else:
            x1 = xO1 - root*z1
            x2 = xO2 - root*z2
            x3 = xO3 - root*z3
            print ("\n")
            print ("x1: ",float(x1))
            print ("x2: ",float(x2))
            print ("x3: ",float(x3))
            return float(x1),float(x2),float(x3)


    def finalGx(self,x1,x2,x3):
        if x3 ==None:
            x = self.f1(x1,x2)
            print ("Final g val: " + str(x[2]))
            return x[2]
        else:
            x = self.f1Version2(x1,x2,x3)
            print ("\nFINAL g VALUE: " + str(x[3]))
            return x[3]

    def initialValueSet(self,a,b,c):
        self.itZero = [a,b,c]

        print ("\n" * 2 + "x(" + str(self.itNumber) + ")" + "\n")
        print (self.itZero)


    def mainLoop(self):
        TolCheck = 999999999999

        while TolCheck > self.tolerance:
            print ("ITERATION NUMBER: " + str(self.itNumber))
            print ("*" * 100)
            qq = self.theJacobian(self.itZero[0],self.itZero[1],self.itZero[2])
            #print "JACOBIAN MATRIX WITHOUT SIMPLIFICATION"
            #print matrix(qq)
            print ("DIAGONALIZED JACOBIAN")
            print (numpy.matrix(numpy.diag(numpy.diag(qq))))
            #print "Echelon form: " + str(qq.echelon_form())

            x = self.zFind(numpy.matrix(numpy.diag(numpy.diag(qq))), self.itZero[0], self.itZero[1], self.itZero[2])
            if x[2] == None:
                k = self.g3(self.alpha3, self.tolerance, numpy.matrix(2,1,[x[0], x[1]]), x[3], self.itZero)
            else:

                k = self.g3(self.alpha3, self.tolerance, numpy.matrix(3,1,[x[0], x[1], x[2]]), x[3], self.itZero)
            z = self.findHVal(k[0],k[1], k[2], k[3], k[4])
            q = self.createNewtonIntPoly(z[0], z[1], z[2], z[3], z[4])
            p = self.findFirstIte(q, self.itZero[0], self.itZero[1], self.itZero[2], x[0], x[1], x[2])

            #Class variable to track iteration numbers
            self.itNumber += 1

            finalG = self.finalGx(p[0],p[1],p[2])
            #TolCheck = finalG^2 - k[0]^2
            TolCheck = abs(finalG - k[0])
            #print "|g^k - g|: " + str(abs(sqrt(TolCheck)))
            print ("|g^k - g|: " + str(abs(TolCheck)))

            #This method assigns the new values to class variable self.itZero
            self.initialValueSet(p[0], p[1], p[2])

            print ("*" * 100)
            print ("\n"*5)
            print ("Exact Answer: " + str([5,0,-0.5235988]))


hhh = SteepestDescents([0,0,0])
hhh.mainLoop()

if __name__ == '__main__':
    pass