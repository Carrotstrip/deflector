from visual import *
import random

"""All units are SI base units."""

rScale = 4
mSun = 1.989e30
rSun = 695.508e6

muEarth = 3.986e14
aEarth = 149.6e9
rPEarth = 147.1e9
eEarth = 0.0167086
rEarth = 6371e3
windowRange = rPEarth*1.5

class Body(sphere):

    GRAVC = 6.674e-11

    def __init__(self, mass, trueRadius, visualRadius):
        sphere.__init__(self, radius = visualRadius)
        self.mass = mass
        self.mu = GRAVC*mass
        self.trueRadius = trueRadius
        #self.velocity = vector(0, 0, 0)
        #self.acceleration = vector(0, 0, 0)
        

class Planet(Body):

    rPlanet = 7e8

    def __init__(self, mass, trueRadius):
        Body.__init__(self, mass, trueRadius, self.rPlanet)
        

class smallBody(Body):

    rSmallBody = 4e8

    def __init__(self, mass, trueRadius):
        Body.__init__(self, mass, trueRadius, self.rSmallBody)
        

class SolarSystem:

    dt = 1

    def __init__(self, star):
        self.star = star
        self.planets = []
        self.smallBodies = []

    def addBody(self, body, velocity, position):
        """Add a body to the SolarSystem. The only body
        that doesn't get added this way is the star."""
        # we define position and velocity here because they
        # have no meaning until put in the context of a SolarSystem
        body.velocity = velocity
        body.position = position
        # we can now, given velocity and position relative
        # to the star, define all orbital parameters, which we will do
        
        # there are only Planets and smallBodies
        if isinstance(body, Planet):
            self.planets.append(body)
        else:
            self.smallBodies.append(body)

    def inSOI(self, body2, body1):
        return

    def getAccelerations(self):
        """Get accelerations for all the objects. Use sphere of influence
        physics for patched conic approximation. This code only allows for
        a body to be in one SOI at a time. A larger body cannot be in a
        smaller body's SOI. Small bodies are bodies that can be inside a
        planet's SOI; planets are bodies which are large enough
        to have an SOI but are always only within the SOI of the star.
        The star is motionless in our frame."""
        for body2 in self.smallBodies:
            body2.acc = -(body2.pos-star.pos)*star.mu/(mag(body2.pos-star.pos)**3)
            for body1 in self.largeBodies:
                if self.inSOI(body2, body1):
                    body2.acc = -(body2.pos-body1.pos)*body1.mu/(mag(body2.pos-body1.pos)**3)

    def updateVelocities(self):
        """Update velocities for all celestial bodies in the SolarSystem."""
        for body in self.smallBodies+self.planets:
            body.velocity += body.acceleration*self.dt

    def updatePositions(self):
        """Update positions for all celestial bodies in the SolarSystem."""
        for body in self.smallBodies+self.planets:
            body.pos += body.velocity*self.dt





            
        
scene = display(title = "Asteroid Deflection", width=800, height=800, range=rPEarth*1.5)
sun = Body(mSun, rSun, rSun*9)
sun.color = color.yellow
home = SolarSystem(sun)
earth = Planet(400, 6371e3, vector(rPEarth, 0, 0))
home.addBody(earth)

while 1:
    rate(10000)
    home.getAccelerations()
    home.updateVelocities()
    home.updatePositions()























