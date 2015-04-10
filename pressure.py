'''Based on http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/GDC03.pdf
/**
 * Copyright (c) 2009 Oliver Hunt <http://nerget.com>
 * 
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */'''
# Note: Check to see ++currentrow doesnt throw errors.
# Note: Check for loops of lines 28 - 190 so equations that use i,j,etc to see if they are used.
# Note: Replaced all (this.) with (self.). Who knows if it will convert over.
import math
class FluidField():
    #create innit thingy similar to physics hw assignments.
    def __init__(self,pos, vel, mass = 1.0, size = 10, shape = "circle", color = (0,0,255)):
        self.iterations = 10 
        self.visc = 0.5 
        self.dt = 0.1 
        self.dens 
        self.dens_prev 
        self.u 
        self.u_prev 
        self.v 
        self.v_prev 
        self.width 
        self.height 
        self.rowSize 
        self.size
        self.displayFunc  
    def addFields(self,x, s, dt):
        i = 0
        for i in range(0,self.size):
            x[i] += dt*s[i]
    

    def set_bnd(self,b, x):
        if b==1:
            
            for i in range(0,self.size): 
                x[i] =  x[i + self.rowSize] 
                x[i + (self.height+1) *self.rowSize] = x[i + self.height * self.rowSize] 
            

            for j in range(0,self.height): 
                x[j * self.rowSize] = -x[1 + j * self.rowSize] 
                x[(self.width + 1) + j * self.rowSize] = -x[self.width + j * self.rowSize]
            
        elif (b == 2):
            for i in range(1,self.width+1):
                x[i] = -x[i + self.rowSize]
                x[i + (self.height + 1) * self.rowSize] = -x[i + self.height * self.rowSize] 
            

            for j in range(1,self.height+1):
                x[j * self.rowSize] =  x[1 + j * self.rowSize] 
                x[(self.width + 1) + j * self.rowSize] =  x[self.width + j * self.rowSize] 
            
        else: 
            for i in range(1, self.width+1):
                x[i] =  x[i + self.rowSize] 
                x[i + (self.height + 1) * self.rowSize] = x[i + self.height * self.rowSize] 
            

            for j in range(1, self.height+1):
                x[j * self.rowSize] =  x[1 + j * self.rowSize] 
                x[(self.width + 1) + j * self.rowSize] =  x[self.width + j * self.rowSize] 
            
        
        maxEdge = (self.height + 1) * self.rowSize 
        x[0]                 = 0.5 * (x[1] + x[self.rowSize]) 
        x[maxEdge]           = 0.5 * (x[1 + maxEdge] + x[self.height * self.rowSize]) 
        x[(self.width+1)]         = 0.5 * (x[self.width] + x[(self.width + 1) + self.rowSize]) 
        x[(self.width+1)+maxEdge] = 0.5 * (x[self.width + maxEdge] + x[(self.width + 1) + self.height * self.rowSize]) 
    

    def lin_solve(self,b, x, x0, a, c):
    
        if a == 0 and c == 1: 
            for j in range(1, self.height+1):
                currentRow = j * self.rowSize 
                ++currentRow 
                for i in range(0,self.width):
                    x[currentRow] = x0[currentRow] 
                    ++currentRow 
                
            
            self.set_bnd(b, x) 
        else:
            invC = 1 / c 
            for k in range(0, self.iterations):
                for j in range(0,self.height):
                    lastRow = (j - 1) * self.rowSize 
                    currentRow = j * self.rowSize 
                    nextRow = (j + 1) * self.rowSize 
                    lastX = x[currentRow] 
                    ++currentRow 
                    for i in range(0,self.width):
                        lastX = x[currentRow] = (x0[currentRow] + a*(lastX+x[++currentRow]+x[++lastRow]+x[++nextRow])) * invC 
                
                self.set_bnd(b, x) 
            
        
    
    
    def diffuse(self,b, x, x0, dt):
    
        a = 0 
        self.lin_solve(b, x, x0, a, 1 + 4*a) 
    
    
    def lin_solve2(self,x, x0, y, y0, a, c):
    
        if (a == 0 and c == 1):
            for j in range(1, self.height+1):
                currentRow = j * self.rowSize 
                ++currentRow 
                for i in range(0,self.width):
                    x[currentRow] = x0[currentRow] 
                    y[currentRow] = y0[currentRow] 
                    ++currentRow 
                
            
            self.set_bnd(1, x) 
            self.set_bnd(2, y) 
        else:
            invC = 1/c 
            for k in range(0, self.iterations):
                for j in range(0, self.height):
                    lastRow = (j - 1) * self.rowSize 
                    currentRow = j * self.rowSize 
                    nextRow = (j + 1) * self.rowSize 
                    lastX = x[currentRow] 
                    lastY = y[currentRow] 
                    ++currentRow 
                    for i in range(0, self.width):
                        lastX = x[currentRow] = (x0[currentRow] + a * (lastX + x[currentRow] + x[lastRow] + x[nextRow])) * invC 
                        lastY = y[currentRow] = (y0[currentRow] + a * (lastY + y[++currentRow] + y[++lastRow] + y[++nextRow])) * invC 
                    
                
                self.set_bnd(1, x) 
                self.set_bnd(2, y) 
            
        
    
    
    def diffuse2(self,x, x0, y, y0, dt):
    
        a = 0 
        self.lin_solve2(x, x0, y, y0, a, 1 + 4 * a) 
    
    
    def advect(self,b, d, d0, u, v, dt):
    
        Wdt0 = dt * self.width 
        Hdt0 = dt * self.height 
        Wp5 = self.width + 0.5 
        Hp5 = self.height + 0.5 
        for j in range(1, self.height+1):
            pos = j * self.rowSize 
            for i in range(1,self.width+1):
                x = i - Wdt0 * u[++pos]  
                y = j - Hdt0 * v[pos] 
                if (x < 0.5):
                    x = 0.5 
                elif (x > Wp5):
                    x = Wp5 
                i0 = x | 0 
                i1 = i0 + 1 
                if (y < 0.5):
                    y = 0.5 
                elif (y > Hp5):
                    y = Hp5 
                j0 = y | 0 
                j1 = j0 + 1 
                s1 = x - i0 
                s0 = 1 - s1 
                t1 = y - j0 
                t0 = 1 - t1 
                row1 = j0 * self.rowSize 
                row2 = j1 * self.rowSize 
                d[pos] = s0 * (t0 * d0[i0 + row1] + t1 * d0[i0 + row2]) + s1 * (t0 * d0[i1 + row1] + t1 * d0[i1 + row2]) 
            
        
        self.set_bnd(b, d) 
    
    
    def project(self,u, v, p, div):
    
        h = -0.5 / math.sqrt(self.width * self.height) 
        for j in range(1,self.height+1):
            row = j * self.rowSize 
            previousRow = (j - 1) * self.rowSize 
            prevValue = row - 1 
            currentRow = row 
            nextValue = row + 1 
            nextRow = (j + 1) * self.rowSize 
            for i in range(1,self.width+1):
                div[++currentRow] = h * (u[++nextValue] - u[++prevValue] + v[++nextRow] - v[++previousRow]) 
                p[currentRow] = 0 
            
        
        self.set_bnd(0, div) 
        self.set_bnd(0, p) 
        
        self.lin_solve(0, p, div, 1, 4 ) 
        wScale = 0.5 * self.width 
        hScale = 0.5 * self.height 
        for j in  range(1, self.height+1):
            prevPos = j * self.rowSize - 1 
            currentPos = j * self.rowSize 
            nextPos = j * self.rowSize + 1 
            prevRow = (j - 1) * self.rowSize 
            currentRow = j * self.rowSize 
            nextRow = (j + 1) * self.rowSize 

            for i in range(1,self.width+1):
                u[++currentPos] -= wScale * (p[++nextPos] - p[++prevPos]) 
                v[currentPos]   -= hScale * (p[++nextRow] - p[++prevRow]) 
            
        
        self.set_bnd(1, u) 
        self.set_bnd(2, v) 
    
    
    def dens_step(self,x, x0, u, v, dt):
    
        self.addFields(x, x0, dt) 
        self.diffuse(0, x0, x, dt ) 
        self.advect(0, x, x0, u, v, dt ) 
    
    
    def vel_step(self,u, v, u0, v0, dt):
    
        self.addFields(u, u0, dt ) 
        self.addFields(v, v0, dt ) 
        temp = u0 
        u0 = u 
        u = temp 
        temp2 = v0  
        v0 = v  
        v = temp2 
        self.diffuse2(u,u0,v,v0, dt) 
        self.project(u, v, u0, v0) 
        temp3 = u0 
        u0 = u 
        u = temp3  
        temp4 = v0 
        v0 = v 
        v = temp4 
        self.advect(1, u, u0, u0, v0, dt) 
        self.advect(2, v, v0, u0, v0, dt) 
        self.project(u, v, u0, v0 ) 
    #confusing bit here. Definitely double check    
    def tempCallback(self,d,u,v): {} 
    #will most likely need to put uiCallback and any other lambda variable in innit so that it is avaialable everywhere. Would fix a few translation errors i am noticing
    uiCallback = lambda(self,d,u,v):self.tempCallback(d,u,v)
 
    def setVelocityTemp(self,x, y, xv, yv): 
        self.u[(x + 1) + (y + 1) * self.rowSize] = xv 
        self.v[(x + 1) + (y + 1) * self.rowSize] = yv 
    def setDensityTemp(self,x,y,d):
        self.dens[(x + 1) + (y + 1) * self.rowSize] = d
    
    def Field(self,dens, u, v):
        #Just exposing the fields here rather than using accessors is a measurable win during display (maybe 5%)
        #but makes the code ugly.
        self.setDensity = lambda(self,x, y, d): self.setDensityTemp(x,y,d)
        
        self.getDensity = lambda(self,x, y):  dens[(x + 1) + (y + 1) * self.rowSize]
        
        self.setVelocity = lambda(self,x, y, xv, yv): self.setVelocityTemp(x,y,xv,yv)
        
        self.getXVelocity = lambda(self,x, y): u[(x + 1) + (y + 1) * self.rowSize] 
        
        self.getYVelocity = lambda(self,x, y): v[(x + 1) + (y + 1) * self.rowSize] 
        
        #Setter and getter using lambda. Double check the whole self.width in the body of lambda.
        self.width = lambda self:  self.width 
        self.height = lambda self:  self.height
    
    def queryUI(self,d, u, v):
    
        for i in range(0, self.size):
            u[i] = v[i] = d[i] = 0.0
        # Supposed to create new instance of object. Python seems to do this auto without new. Double check.
        # Causes error if new is before Field    
        self.uiCallback(self.Field(d, u, v))
    
    def updateDef(self):
        self.queryUI(self.dens_prev, self.u_prev, self.v_prev) 
        self.vel_step(self.u, self.v, self.u_prev, self.v_prev, self.dt) 
        self.dens_step(self.dens, self.dens_prev, self.u, self.v, self.dt) 
        self.displayFunc(self.Field(self.dens, self.u, self.v)) 
    
    #original code for line 313:  this.update = function () { queryUI(dens_prev, u_prev, v_prev)  and so on   
    update = lambda (self):self.updateDef()
    
    def setDisplayFunctionTemp(self,func):
        self.displayFunc = func
    setDisplayFunction = lambda(self,func): self.setDisplayFunctionTemp(func)
    
    #another getter and possible innit variable
    iterations = lambda(self): { self.iterations }
    
    def setIterationz(self,iters):
        if (iters > 0 and iters <= 100):
           self.iterations = iters 
    setIterations = lambda(self,iters): self.setIterationz(iters)
    
    # Not sure if this is neccessary. Though looks like this is where they set up the frame in the html file.
    def setuiCallBackplz(self,callback):
        uiCallback = callback
        
    #possible innit variable    
    setUICallback = lambda(self,callback):self.setuiCallBackplz(callback)
   
    def reset(self):
    
        self.rowSize = self.width + 2 
        self.size = (self.width+2)*(self.height+2) 
        self.dens = []*self.size 
        self.dens_prev = []*self.size 
        self.u = []*self.size 
        self.u_prev = []*self.size 
        self.v = []*self.size 
        self.v_prev = []*self.size 
        for i in range(0,self.size):
            self.dens_prev[i] = self.u_prev[i] = self.v_prev[i] = self.dens[i] = self.u[i] = self.v[i] = 0 
            
    # The line below will need to be changed or more likely moved to a different function to get triggered
    #reset = reset()
    
    def setNewResolution(self,hRes, wRes):
        res = wRes * hRes 
        if (res > 0 and res < 1000000 and (wRes != self.width or hRes != self.height)):
            self.width = wRes 
            self.height = hRes 
            self.reset() 
            return True
        
        return False 
        
    #possible innit variables
    setResolution = lambda(self,hRes,wRes): self.setNewResoluation(self,hRes,wRes)
    
    #problem with this line. says need only 1 arguement, when it should need 2.
    #setResolution(64,64) 


print "win"