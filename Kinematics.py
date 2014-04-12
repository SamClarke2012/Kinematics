import math
import time 
## @package Kinematics
#  A rigid body kinematics library for python. Holds information about the attitude
#  of a rigid body relative to a target attitude.

##
#  A Quaternion object with built-in functions to convert to Euler angles (Tait-Bryan) 
class Quaternion():

    ##
    #  A Quaternion object with built-in functions to work with quaternions in order to 
    #  add rotations, find the rotational difference between two quaternions, and conversion
    #  to quaternions to Euler angles(phi,theta,psi)(z,y,x)(yaw, pitch, roll).
    #
    # @param self The object pointer
    # @param w The \f$W\f$ component of \f$q = W+Xi+Yj+Zk\f$
    # @param x The \f$X\f$ component of \f$q = W+Xi+Yj+Zk\f$
    # @param y The \f$Y\f$ component of \f$q = W+Xi+Yj+Zk\f$
    # @param z The \f$Z\f$ component of \f$q = W+Xi+Yj+Zk\f$
    def __init__(self, w=1, x=0, y=0, z=0):
        a = [w,x,y,z]
        for i in a:
            if isinstance(i, int):
                i = float(i)
        self.w = w
        self.x = x
        self.y = y
        self.z = z
        self.normalise()
 
    def __str__(self):
        return str(self.w)+'+'+str(self.x)+'i+'+str(self.y)+'j+'+str(self.z)+'k'
 
    ##
    #  Convert the current attitude to Euler angles (Tait-Bryan ZYX)\n
    #  Yaw   = (\f$\phi\f$ or rotation arouns the \f$Z\f$ axis )  = \f$\arctan(2.0(Q(w)Q(z)+Q(x)Q(y)),1.0-2.0(Q(y)^2+Q(z)^2))\f$\n
    #  Pitch = (\f$\theta\f$ or rotation arouns the \f$Y\f$ axis )  = \f$\arcsin(2.0(Q(w)Q(y)-Q(z)Q(x)))\f$\n
    #  Roll  = (\f$\psi\f$ or rotation arouns the \f$X\f$ axis )  = \f$\arctan(2.0(Q(w)Q(z)+Q(x)Q(y)),1.0-2.0(Q(y)^2+Q(z)^2))\f$\n
    #  @param Returns A EulerAngles object
    def toEuler(self):
        """
        Returns a EulerAngles object of itself (in degress)
        
        """                
        # yaw = atan2(2.0f*(qw*qz + qx*qy), 1-2.0*(qy*qy + qz*qz))
        # pitch = wrap_asin(2.0*(qw*qy - qz*qx))
        # roll = (atan2(2.0*(qw*qx + qy*qz), 1.0 - 2.0*(qx*qx + qy*qy)))
        yaw = math.atan2(2.0*(self.w*self.z + self.x*self.y),\
                        (1.0 - 2.0*(self.y**2 + self.z**2)))/math.pi*180.0

        roll = math.atan2(2.0*(self.w*self.x + self.y*self.z), \
                         1.0 - 2.0*(self.x**2 + self.y**2))/math.pi*180.0

        pitch = _wrap_asin(2.0*(self.w*self.y - self.z*self.x))/math.pi*180.0

        return EulerAngles(yaw, pitch, roll)
 
    ##
    # Find the product of this quaternion and another quaternion.
    # @param q A Quaternion object that will be used to solve to product
    # @param Returns A new quaternion as the result
    def magnitude(self):
        """
        Returns the magnitude of this quaternion
        """
        return math.sqrt(self.w**2 + self.x**2 + self.y**2 + self.z**2)
 
    ##
    # Normalise the components of this quaternion.
    def normalise(self):
        """
        Normalise this quaternion
        """
        try:
            m = self.magnitude()
            self.w /= m
            self.x /= m
            self.y /= m
            self.z /= m
        except ZeroDivisionError:
            # Quaternion was init with no values
            pass

    ##
    # Find the conjugate of this quaternion
    # @param Returns A new quaternion as the result 
    def conjugate(self):
        """
        Returns a new quaternion that is the conjugate of this quaternion
        """
        return Quaternion(self.w, -self.x, -self.y, -self.z)
 
    ##
    # Find the quaternion difference between this quaternion and another.
    # @param q A Quaternion object that will be used
    # @param Returns A new quaternion as the result\n
    def difference(self, q):
        """
        q = Quaternion
        Returns a new quaternion that is difference between this quaternion and q
        """
        # DeltaQ = Q(dest) * Q(origin)^{-1}
        return q.productOf(self.conjugate())
    ##
    # Find the product of this quaternion and another quaternion.
    # @param q A Quaternion object that will be used to solve to product
    # @param Returns A new quaternion as the result
    def productOf(self, q):
        """
        q = Quaternion
        returns a quaternion thst is the product of this quaternion and q
        """               # (Q1 * Q2).w = (w1w2 - x1x2 - y1y2 - z1z2)
        return Quaternion(self.w*q.w - self.x*q.x - self.y*q.y - self.z*q.z,\
                          # (Q1 * Q2).x = (w1x2 + x1w2 + y1z2 - z1y2)
                          self.w*q.x + self.x*q.w + self.y*q.z - self.z*q.y,\
                          # (Q1 * Q2).y = (w1y2 - x1z2 + y1w2 + z1x2)
                          self.w*q.y - self.x*q.z + self.y*q.w + self.z*q.x,\
                          # (Q1 * Q2).z = (w1z2 + x1y2 - y1x2 + z1w2)
                          self.w*q.z + self.x*q.y - self.y*q.x + self.z*q.w)


##
#  A Euler angles object with built-in functions to convert to quaternions 
class EulerAngles():
    """
        A basic Euler object
    """
    def __init__(self, yaw = 0, pitch = 0, roll = 0):
        for i in [yaw, pitch, roll]:
            if isinstance(i, int):
                i = float(i)
        self.yaw = yaw
        self.pitch = pitch
        self.roll = roll      
 
    def __str__(self):
        return 'Yaw(z): '+str(self.yaw)+' Pitch(y): '\
                +str(self.pitch)+' Roll(x): '+str(self.roll)
 
    def getYaw(self):
        return self.yaw
 
    def getPitch(self):
        return self.pitch
 
    def getRoll(self):
        return self.roll
    ##
    # Return a vector of this attitude
    # @param Returns A tuple (Y,P,R) as the result\n 
    def YPRvector(self):
        return (self.yaw, self.pitch, self.roll)
    ##
    # Convert to a quaternion
    # @param Returns A new quaternion as the result\n
    def toQuaternion(self):
        """
        Returns a quaternion object of itself
        """  
        cosRoll  = math.cos(self.roll*math.pi/180.0*0.5)
        sinRoll  = math.sin(self.roll*math.pi/180.0*0.5)
        cosPitch = math.cos(self.pitch*math.pi/180.0*0.5)
        sinPitch = math.sin(self.pitch*math.pi/180.0*0.5)
        cosYaw   = math.cos(self.yaw*math.pi/180.0*0.5)
        sinYaw   = math.sin(self.yaw*math.pi/180.0*0.5)

        w = cosRoll*cosPitch*cosYaw + sinRoll*sinPitch*sinYaw
        x = sinRoll*cosPitch*cosYaw - cosRoll*sinPitch*sinYaw
        y = cosRoll*sinPitch*cosYaw + sinRoll*cosPitch*sinYaw
        z = cosRoll*cosPitch*sinYaw - sinRoll*sinPitch*cosYaw

        return Quaternion(w,x,y,z)

 
##
#  A Kinematics object that contains the attitude of the rigid body,\n
#  a target attitude, the angular velocity of the rigid body, and the error\n
#  of the body (target - self) 
class Kinematics(Quaternion):
    def __init__(self,  w=1, x=0, y=0, z=0):
        Quaternion.__init__(self, w, x, y, z)
        self.target = self
        self.time = time.time()
        self.degSecYPR = (0.0, 0.0, 0.0)
        self.error = self.difference(self.target)
 
    def __str__(self):
        data = 'Current Quaternion (W+Xi+Yj+Zk):\n'
        data += str(self.w)+'+'+str(self.x)+'i+'+str(self.y)+'j+'+str(self.z)+'k\n\n'
        data += 'Current attitude:\n'
        data += str(self.toEuler())+'\n\n'
        data += 'Target attitude:\n'
        data += str(self.target.toEuler())+'\n\n'
        data += 'Error:\n'
        data += str(self.error.toEuler())+'\n\n'
        data += 'Angular velocities (Deg/Sec):\n'
        data += str(self.getAngularVelocities())+'\n'
        return data
    ##
    # Set an new target quaternion
    # @param q A Quaternion to use as the target attitude
    def setTarget(self, q):
        """
        Set new target attitude
        """
        self.target = q
        self.update(self, q)
    ##
    # Rotate current target by a quaternion
    # @param q A Quaternion to use as the modifier
    def rotateTarget(self, q):
        """
        Rotate the target attitude by a quaternion
        """
        self.target = self.target.productOf(q)
        self.update(self, self.target)
    ##
    # Get current angular velocities as a vector (Y,P,R)
    # @param Returns A tuple containing the velocity (deg/sec) for each axis (Y,P,R)
    def getAngularVelocities(self):
        """
        Get angular velocity in (Y,P,R)deg/sec
        """
        return self.degSecYPR
    ##
    # Get current self/target rotational error
    # @param Returns A Quaternion object
    def getErrorQuaternion(self):
        """
        Get error (self vs target) as a quaternion
        """
        return self.error
    ##
    # Get current self/target rotational error
    # @param Returns A EulerAngles object
    def getErrorEuler(self):
        """
        Get error (self vs target) as a Euler object
        """
        return self.error.toEuler()
    ##
    # Update current attitude and target attitude, calculate angular velocities\n
    # and realtime error in attitude.
    # @param attitude A Quaternion object representing our current attitude
    # @param target A Quaternion object representing a target attitude
    def update(self, attitude, target):
        """
        Update current attitude (attitude) and target attitude (target),
        calculate angular velocities and realtime error.
        """
        # delta T
        prevTime = self.time
        self.time = time.time()
        dt = self.time - prevTime
        # store attitude
        prev = Quaternion(self.w, self.x, self.y, self.z)
        # adopt new attitude
        self.w = attitude.w
        self.x = attitude.x
        self.y = attitude.y
        self.z = attitude.z
        # adopt new target rotation
        self.target = target
        # calculate attitude change as Euler angles
        delta = prev.difference(self).toEuler()
        # calculate angular velocities in deg/sec
        self.degSecYPR = (delta.getYaw()/dt, delta.getPitch()/dt, delta.getRoll()/dt)
        # calculate target error as a unit quaternion
        self.error = self.difference(self.target)
 
       
#####################################################################################
 
         
def _wrap_asin(v):
    """
    Constrain asin inputs
    """
    if v >= 1.0: return math.pi/2.0
    elif v <= -1.0: return -math.pi/2.0
    else: return math.asin(v)
    

#####################################################################################    
    
    #  Tests

#####################################################################################
##    
##if __name__=="__main__":
##    
##    def test1():
##        # Euler to quaternion back to Euler
##        for n in xrange(91):
##            # Start in good ol' Euler angles
##            e = EulerAngles(n*2,n-1,n/2.0)
##            q = e.toQuaternion()
##            print '#'*3+'TEST CASE = ',e,'#'*3
##            print 'To Quaternion = ', q
##            print 'Back to Euler = ', q.toEuler(), '\n'
##     
##    def test2():
##        # Target error test
##        k = Kinematics() 
##        # Yaw 0, Pitch 0, Roll +90
##        k.update(k, Quaternion(1,1,0,0))
##        print k
##      
##    def test3():
##        # Make a relative modification to the target attitude  
##        # (relative change described as both Euler and Quaternion)
##        body = Kinematics()
##        print body
##        print '#'*30,'\nAdding +89 degrees pitch...\n','#'*30
##        # As a conversion from Euler angles
##        change = EulerAngles(0,89,0).toQuaternion()
##        body.rotateTarget(change)
##        print body
##        print '#'*30,'\nAdding +90 degrees roll...\n','#'*30
##        # As a quaternion
##        change = Quaternion(1,1,0,0)
##        body.rotateTarget(change)
##        print body
        
