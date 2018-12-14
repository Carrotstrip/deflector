from vpython import *
import numpy as np
import random
import orbital
import time
from poliastro import iod
from astropy import units as u

"""This code generates an asteroid of random mass and radius and solves the Lambert
Problem to give it the proper initial velocity to intercept Earth in exactly advanceNotice years.
It then generates a kinetic impactor to intercept the asteroid and impart a delta-v,
thus saving the Earth. The following assumptions are made:

Humanity discovers the fatal trajectory on the vernal equinox, aka periapsis.

Lucky for humanity, they had their defenses all set and can launch them on day zero.

The asteroid's orbital plane is the same as Earth's.

"""

"""Thoughts for me:

Can I control where the kinetic impactor hits the asteroid?

If I can't, does it matter? The constraint here is TOF in Lambert,
I'd have to set that to the proper value that corresponds with where the probe hits 
the asteroid, but how do I find that value?

I should try Earth-based lasers instead that'd be easier. I could just impart a constant acceleration
on the asteroid, bing bang boom done. How to make sure it's enough?

How do I make sure either nudge is enough? Can I? With the kinetic impactor, I know
vf before the sim runs, so I can iterate on different probe masses in that direction until
I find one that will result in non collision.

What is the non collision condition without running the sim though? How do I analytically get perturbation?
Both orbits are fully defined, I should be able to right?

Idea to get around that: have the user put in all the relevant values (mass, radius, advanceNotice, mass of
kinetic impactor, laser thrust imparted) and then just try and see if we can deflect it.

I could improve the kinetic impactor by making the collision as head on as possible, how?

Now we should make a database of a few IRL asteroids and see what it takes to deflect them given certain
advance notices. Use both lasers and kinetic impactors

We should also shift the generateAsteroid functionality into the more general generateInterceptor
generateAsteroid only works to intercept bodies that are at periapsis at the epoch; I will leave
it intact somewhere in case anybody wants to see what I mean

"""

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

    def __init__(self, name, mass, trueRadius, visualRadius):
        sphere.__init__(self, radius = visualRadius, make_trail = False)
        self.name = name
        self.color = color.red
        self.mass = mass
        self.mu = self.GRAVC*mass
        self.trueRadius = trueRadius
        self.acceleration = vector(0, 0, 0)
        self.velocity = vector(0, 0, 0)
        

class Planet(Body):

    rPlanet = 30e8

    def __init__(self, name, mass, trueRadius):
        Body.__init__(self, name, mass, trueRadius, self.rPlanet)
        

class smallBody(Body):

    rSmallBody = 20e8

    def __init__(self, name, mass, trueRadius):
        Body.__init__(self, name, mass, trueRadius, self.rSmallBody)
        

class SolarSystem:

    dt = 100
    t = 0
    velArr = False

    def __init__(self, star):
        self.star = star
        self.planets = {}
        self.smallBodies = {}

    def addBody(self, body, velocity, position):
        """Add a body to the SolarSystem. The only body
        that doesn't get added this way is the star."""
        # we define position and velocity here because they
        # have no meaning until put in the context of a SolarSystem
        body.velocity = vector(velocity[0], velocity[1], velocity[2])
        body.pos = vector(position[0], position[1], position[2])
        # we can now, given velocity and position relative
        # to the star, define all orbital parameters, which we will do
        position = np.array(position, dtype=np.float32)
        velocity = np.array(velocity, dtype=np.float32)
        body.elements = orbital.utilities.elements_from_state_vector(position, velocity, self.star.mu)
        body.h = orbital.angular_momentum(position, velocity)
        # there are only Planets and smallBodies
        if isinstance(body, Planet):
            self.planets[body.name] = body
        else:
            self.smallBodies[body.name] = body

    def getBody():
        pass

    # def generateAsteroid(self, name, planetToHit):
    #     """Generates an asteroid on a collision course with Earth.
    #     this requires that the orbits are at the same place at the same time."""
    #     #mass = raw_input('Asteroid Mass: ')
    #     #advanceNotice = raw_input('Time to Impact: ')
    #     mass = 78e9#random.randint(2e14, 6e15)
    #     advanceNotice = .3#random.randint(5, 10)
    #     # random radius
    #     trueRadius = 5000#random.randint(5, 15)
    #     # velocity and position have to be such that the asteroid collides with Earth
    #     # Fix collision location?
    #     # start with a position then we solve the 'given location, find time till it gets there' problem
    #     # and then give it the proper velocity such that it gets there when earth does
    #     # take time till impact and find out where earth will be when that happens, given
    #     # that it starts from periapsis
    #     # then take the asteroid's starting position and give it the requisite velocity to
    #     # intersect the collision point at tImpact
    #     # you know two points in an orbit and the time it takes to get between them, you can then figure out
    #     # the velocity vector at rStart
    #     tMinus = yearsToSeconds(advanceNotice)
    #     # M is how many radians of mean motion till impact
    #     M = sqrt(self.star.mu/(planetToHit.elements.a)**3)*(tMinus)
    #     # theta is true anomaly of planet at tImpact
    #     theta = orbital.utilities.true_anomaly_from_mean(planetToHit.elements.e, M)
    #     r = (np.linalg.norm(planetToHit.h)**2/self.star.mu)/(1+planetToHit.elements.e*cos(theta))
    #     rImpact = [r*cos(theta), r*sin(theta), 0]
    #     # s = sphere(pos=vector(rImpact[0], rImpact[1], rImpact[2]), radius = rSun*6)
    #     asteroid = smallBody(name, mass, trueRadius)
    #     r = [rPEarth, rPEarth, 0] #[random.randint(-rPEarth*1.1, rPEarth*1.1), random.randint(rPEarth*1.1, rPEarth*1.3), 0]
    #     # now what is the v that takes it through rImpact at tImpact?
    #     # this is lambert's problem!
    #     (v0, v) = self.lambert(r, rImpact, tMinus)
    #     self.addBody(asteroid, v0, r)

    def generateInterceptor(self, body, bodyToIntercept, initPosition, tImpact):
        # mass = body.mass
        # trueRadius = body.radius
        tMinus = yearsToSeconds(tImpact)
        # TODO this doesn't work since the asteroid doesn't necessarily start from its periapsis
        # i need to add the mean anomaly of the asteroid at t = 0 (sim start), how do i get this
        # use the f
        true = bodyToIntercept.elements.f
        M0 = orbital.mean_anomaly_from_true(bodyToIntercept.elements.e, true)
        M = sqrt(self.star.mu/(bodyToIntercept.elements.a)**3)*(tMinus) + M0
        # theta is true anomaly of bodyToIntercept at tImpact
        theta = orbital.utilities.true_anomaly_from_mean(bodyToIntercept.elements.e, M)
        # print('theta', theta)
        r = (np.linalg.norm(bodyToIntercept.h)**2/self.star.mu)/(1+bodyToIntercept.elements.e*cos(theta))
        # now rotate that into our periapsis on right frame
        position = [bodyToIntercept.pos.x, bodyToIntercept.pos.y, bodyToIntercept.pos.z]
        velocity = [bodyToIntercept.velocity.x, bodyToIntercept.velocity.y, bodyToIntercept.velocity.z]
        position = np.array(position, dtype=np.float32)
        velocity = np.array(velocity, dtype=np.float32)
        right = np.array([1, 0, 0], dtype=np.float32)
        ev = orbital.eccentricity_vector(position, velocity, self.star.mu)
        # a = arrow(pos=vector(0, 0, 0), axis=vector(1e11*ev.x, 1e11*ev.y, 1e11*ev.z), shaftwidth=4e8)
        if ev.y < 0:
            theta -= acos(np.dot(ev, right)/(np.linalg.norm(ev)*np.linalg.norm(right)))
        else:
            theta += acos(np.dot(ev, right)/(np.linalg.norm(ev)*np.linalg.norm(right)))
        rImpact = [r*cos(theta), r*sin(theta), 0]
        # s = sphere(pos=vector(rImpact[0], rImpact[1], rImpact[2]), radius = rSun*6)
        r = [self.planets['earth'].pos.x, self.planets['earth'].pos.y, self.planets['earth'].pos.z]
        # now what is the v that takes it through rImpact at tImpact?
        # this is lambert's problem
        (v0, v) = self.lambert(initPosition, rImpact, tMinus)
        self.addBody(body, v0, initPosition)


    def lambert(self, r0, rf, tof, m=0):
        k = self.star.mu*u.m**3/u.s**2
        r0 = r0*u.m
        rf = rf*u.m
        tof = tof*u.s
        (v0, v), = iod.lambert(k, r0, rf, tof, M=m)
        v0 = [1000*v0.value[0], 1000*v0.value[1], 1000*v0.value[2]]
        v =  [1000*v.value[0], 1000*v.value[1], 1000*v.value[2]]
        return v0, v

    def inSOI(self, smallBody, planet):
        return mag(smallBody.pos - planet.pos) <= planet.elements.a*(planet.mass/self.star.mass)**(2.0/5)

    def fireLasers(self):
        laserThrust = 10
        laserAcc = laserThrust/self.smallBodies['killer'].mass
        laserDirVec = (self.smallBodies['killer'].pos - self.planets['earth'].pos)
        laserAccVec = laserAcc*laserDirVec/mag(laserDirVec)
        # self.laserBeam.pos = self.planets['earth'].pos
        # self.laserBeam.axis = laserDirVec
        # self.laserBeam.shaftwidth = 4e8
        return laserAccVec

    def getAccelerations(self):
        """Get accelerations for all the objects. Use sphere of influence
        physics for patched conic approximation. This code only allows for
        a body to be in one SOI at a time. A larger body cannot be in a
        smaller body's SOI. Small bodies are bodies that can be inside a
        planet's SOI; planets are bodies which are large enough
        to have an SOI but are always only within the SOI of the star.
        The star is motionless in our frame."""
        for body2 in self.smallBodies.values():
            # print(self.fireLasers())
            body2.acceleration = -(body2.pos-self.star.pos)*self.star.mu/(mag(body2.pos-self.star.pos)**3)#+self.fireLasers()
            for body1 in self.planets.values():
                if self.inSOI(body2, body1):
                    self.star.color = color.red
                    # TODO light up the planet that has a smallBody in its SOI
                    #body2.acceleration = -(body2.pos-body1.pos)*body1.mu/(mag(body2.pos-body1.pos)**3)
                else:
                    self.star.color = color.yellow
        for body in self.planets.values():
            body.acceleration = -(body.pos-self.star.pos)*self.star.mu/(mag(body.pos-self.star.pos)**3)

    def updateVelocities(self):
        """Update velocities for all celestial bodies in the SolarSystem."""
        if doCollide(self.smallBodies['killer'], self.smallBodies['savior']):
            # print('boom')
            if not self.velArr:
                self.velArr = True
                inelasticCollision(self.smallBodies['killer'], self.smallBodies['savior'])
        for body in list(self.smallBodies.values())+list(self.planets.values()):
            # print(body.name, body.velocity)
            body.velocity += body.acceleration*self.dt

    def updatePositions(self):
        """Update positions for all celestial bodies in the SolarSystem."""
        for body in list(self.smallBodies.values())+list(self.planets.values()):
            body.pos += body.velocity*self.dt



def doCollide(obj1, obj2):
    """Determine if two objects have collided."""
    # both objects have to be spheres
    return mag(obj1.pos - obj2.pos) <= (obj1.trueRadius+obj2.trueRadius)

def inelasticCollision(obj1, obj2):
    """obj 1 is larger."""
    vec = (obj1.velocity - obj2.velocity)
    obj1.velocity = (obj2.mass*obj2.velocity+obj1.mass*obj1.velocity)/(obj2.mass+obj1.mass)
    obj2.velocity = obj1.velocity

def yearsToSeconds(years):
    return years*525600*60

def vecToList(v):
    return [v.x, v.y, v.z]


scene = canvas(title = "Asteroid Deflection", width=800, height=640, range=windowRange)

sun = Body('sun', mSun, rSun, rSun*9)
sun.color = color.yellow
home = SolarSystem(sun)
earth = Planet('earth', 5.972e24, 6371e3)
earth.color = color.blue
home.addBody(earth, (0, sqrt(home.star.mu/rPEarth), 0), (rPEarth, 0, 0))
killer = smallBody('killer', 488e9, 2e8)
killer.color = color.red
# def generateInterceptor(self, body, bodyToIntercept, initPosition, tImpact)
home.generateInterceptor(killer, earth, [earth.pos.x, earth.pos.x, 0], .4)
savior = smallBody('savior', 1e14, 1e7)
savior.color = color.green
home.generateInterceptor(savior, home.smallBodies['killer'], vecToList(earth.pos), .2)

# scene.bind('keydown', home.generateAsteroid)
years = 0
while not doCollide(home.planets['earth'], home.smallBodies['killer']) and years < 20:
    rate(5000000)  
    home.getAccelerations()
    home.updateVelocities()
    home.updatePositions()
    home.t += home.dt

    # years = float(home.t)/float(525600*60)
if not years >= 20:
    L = label(pos=vector(rPEarth*1.4, rPEarth*1.4, 0),
        text=('Humanity destroyed'), space=30,
        height=16,
        font='sans')
    debris = []
    body = home.planets['earth']
    body.visible = False
    for i in range(20):
        debris.append(sphere(pos=body.pos, radius = body.radius/3))
        vxl = int(body.velocity.x-body.velocity.z/4)
        vxu = int(body.velocity.x*2+body.velocity.z/4)
        vxS = sorted([vxl, vxu])
        vyl = int(body.velocity.y-body.velocity.x/4)
        vyu = int(body.velocity.y*2+body.velocity.x/4)
        vyS = sorted([vyl, vyu])
        vx = random.randint(vxS[0], vxS[1])
        vy = random.randint(vyS[0], vyS[1])
        vz = 0
        debris[-1].velocity = vector(vx, vy, vz)
        debris[-1].color = color.red

    t = 0
    dt = .001
    while t < 100:
        rate(15000)
        for deb in debris:
            deb.pos += deb.velocity*home.dt
        t += dt

"""Now what do we want to know? What useful data can we extract from a run of this sim?

How much did we nudge the asteroid off course at its closest point to Earth?

Can we show a 'ghost' of what would've happened had we not intervened?


"""























