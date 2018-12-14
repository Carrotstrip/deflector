from vpython import *
# from math import *
import random

rScale = 4
muSun = 1.327e20
rSun = 695.508e6

muEarth = 3.986e14
aEarth = 149.6e9
rPEarth = 147.1e9
eEarth = 0.0167086
rEarth = 6371e3
windowRange = rPEarth*1.5

class Body:

    def __init__(self, mass, radius, visualRadius):
        self.mass = mass
        self.radius = radius
        self.visualRadius = visualRadius

class Planet(Body):

    rPlanet = 7e8

    def __init__(self, mass, radius, mu):
        Body.__init__(self, mass, radius, rPlanet)
        self.mu = mu

class smallBody(Body):

    rSmallBody = 4e8

    def __init__(self, mass, radius, visualRadius):
        Body.__init__(self, mass, radius, rSmallBody)
        

class SolarSystem:

    def __init__(self, star):
        self.star = star
        self.largeBodies = []
        self.smallBodies = []

def generateAsteroid():
    """Generates an asteroid on a collision course with Earth.
    this requires that the orbits are at the same place at the same time."""
    global asteroid
    #mass = raw_input('Asteroid Mass: ')
    #advanceNotice = raw_input('Time to Impact: ')
    mass = 500
    advanceNotice = 5
    # random radius
    radius = random.randint(rSun, rSun*5)
    # velocity and position have to be such that the asteroid collides with Earth
    # Fix collision location?
    # start with a position then we solve the 'given location, find time till it gets there' problem
    # and then give it the proper velocity such that it gets there when earth does
    collisionPosition = vector(0, 0, 0)
    position = vector(random.randint(-rPEarth*1.1, rPEarth*1.1), random.randint(-rPEarth*1.1, rPEarth*1.1), 0)
    velocity = vector(-mag(vEarth), 0, 0)
    asteroid.pos = position
    asteroid.radius = radius
    asteroid.mass = mass
    asteroid.velocity = velocity

def doCollide(obj1, obj2):
    """Determine if two objects have collided."""
    # both objects have to be spheres
    return mag(obj1.pos - obj2.pos) <= (obj1.radius+obj2.radius)

def getAccelerations(self):
    """Get accelerations for all the objects. Use sphere of influence
    physics for patched conic approximation. This code only allows for a body to be
    in one SOI at a time. A larger body cannot be in a smaller body's SOI. Small bodies are
    bodies that can be inside a planet's SOI; planets are bodies which are large enough
    to have an SOI but are always only within the SOI of the star. The star is motionless in our frame."""
    for body2 in self.smallBodies:
        body2.acc = -(body2.pos-star.pos)*star.mu/(mag(body2.pos-star.pos)**3)
        for body1 in self.largeBodies:
            if self.inSOI(body2, body1):
                body2.acc = -(body2.pos-body1.pos)*body1.mu/(mag(body2.pos-body1.pos)**3)
        
def updateVelocities(self):
    """Update velocities for all celestial bodies in the sim."""
    for body in self.smallBodies+self.largeBodies:
        body.velocity += body.acceleration*self.dt

def updatePositions(self):
    """Update Positions for all celestial bodies in the sim."""
    for body in self.smallBodies+self.largeBodies:
        body.pos += body.velocity*self.dt

scene = canvas(title = "Asteroid Deflection", width=800, height=800, range=rPEarth*1.5)

# create assets (sun, earth, asteroid, savior)
sun = sphere(color = color.yellow, radius = rSun*rScale*2)
sun.mu = muSun 
earth = sphere(pos = vector(rPEarth,0,0), color = color.green, radius = rSun*rScale)
vEarth = vector(0, sqrt(sun.mu*(2/rPEarth-1/aEarth)), 0)
parkingRad = rEarth + 200e3
savior = sphere(pos = vector(earth.pos.x + parkingRad, earth.pos.y, 0), color = color.cyan, radius = rSun)
savior.velocity = vector(0, sqrt(muEarth/parkingRad), 0)

# make times proportional
tScale = 1500

# set up simulation
dt = 1*tScale
t = 0

# initialize asteroid
#asteroid = sphere(pos = vector(rpEarth, rPEarth, 0), color = color.red, radius = rSun)
#asteroid.mass = 100
#asteroid.velocity = vector(-sEarth*.4, 0, 0)
asteroid = sphere()
generateAsteroid()

#print asteroid

scene.bind('keydown', generateAsteroid)

while not(doCollide(earth, asteroid)):
    rate(1*tScale)
    # get accelerations
    accEarth = -(earth.pos-sun.pos)*sun.mu/(mag(earth.pos-sun.pos)**3)
    accAsteroid = -(asteroid.pos-sun.pos)*sun.mu/(mag(asteroid.pos-sun.pos)**3)
    accSavior = -(savior.pos-earth.pos)*sun.mu/(mag(savior.pos-earth.pos)**3)
    # update velocities
    vEarth += dt*accEarth
    asteroid.velocity += dt*accAsteroid
    savior.velocity = vEarth*1.1# + savior.velocity + dt*accSavior
    # update positions
    earth.pos += dt*vEarth
    asteroid.pos += dt*asteroid.velocity
    savior.pos += dt*savior.velocity
    t += dt

        





        
