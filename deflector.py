from vpython import *
import numpy as np
import random
import orbital
import time
from poliastro import iod
from astropy import units as u
import matplotlib.pyplot as plt

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

class Body(sphere):

    GRAVC = 6.674e-11

    def __init__(self, name, mass, trueRadius, visualRadius):
        sphere.__init__(self, radius = visualRadius, make_trail = True)
        self.name = name
        self.color = color.red
        self.mass = mass
        self.mu = self.GRAVC*mass
        self.trueRadius = trueRadius
        self.acceleration = vector(0, 0, 0)
        self.velocity = vector(0, 0, 0)
        

class Planet(Body):

    rPlanet = 40e8

    def __init__(self, name, mass, trueRadius):
        Body.__init__(self, name, mass, trueRadius, self.rPlanet)
        # find the nearest a smallBody ever came to this planet
        self.nearestHit = 0
        

class smallBody(Body):

    rSmallBody = 30e8

    def __init__(self, name, mass, trueRadius):
        Body.__init__(self, name, mass, trueRadius, self.rSmallBody)
        # did this already undergo its inelastic collision? can only happen once per smallBody
        self.didCollide = False
        # ghosts interact only with gravity
        self.ghost = False
        # set this to true to make this smallBody interact with nonImpulsive forces
        self.nonImpulsiveThrust = 0
        

class SolarSystem:

    dt = 100
    t = 0

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

    def getBody(self, bodyName):
        return ({**(self.smallBodies), **(self.planets)}[bodyName])
    

    def generateInterceptor(self, body, bodyToIntercept, initPosition, tImpact, numRevs):
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
        (v0, v) = self.lambert(initPosition, rImpact, tMinus, numRevs)
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

    def nonImpulsive(self, body, dirVec):
        """Simulates any non impulsive acceleration."""
        # ghosts don't feel this force
        if body.nonImpulsiveThrust <= 0 or body.ghost:
            return vector(0, 0, 0)
        acc = body.nonImpulsiveThrust/body.mass
        accVec = acc*dirVec/mag(dirVec)
        return accVec

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
            body2.acceleration = -(body2.pos-self.star.pos)*self.star.mu/(mag(body2.pos-self.star.pos)**3)+self.nonImpulsive(body2, (body2.pos-self.planets['earth'].pos))
            # for body1 in self.planets.values():
            #     if self.inSOI(body2, body1):
            #         self.star.color = color.red
            #         # TODO light up the planet that has a smallBody in its SOI
            #         #body2.acceleration = -(body2.pos-body1.pos)*body1.mu/(mag(body2.pos-body1.pos)**3)
            #     else:
            #         self.star.color = color.yellow
        for body in self.planets.values():
            body.acceleration = -(body.pos-self.star.pos)*self.star.mu/(mag(body.pos-self.star.pos)**3)

    def updateVelocities(self):
        """Update velocities for all celestial bodies in the SolarSystem."""
        if 'savior' in self.smallBodies:
            if doCollide(self.smallBodies['killer'], self.smallBodies['savior']):
                if not self.smallBodies['savior'].didCollide:
                    # print('boom')
                    self.smallBodies['savior'].didCollide = True
                    # inelasticCollision simulates an inelastic collision between these bodies, 
                    # updating the velocity and sticking them together
                    inelasticCollision(self.smallBodies['killer'], self.smallBodies['savior'])
        for body in list(self.smallBodies.values())+list(self.planets.values()):
            body.velocity += body.acceleration*self.dt

    def updatePositions(self):
        """Update positions for all celestial bodies in the SolarSystem."""
        for body in list(self.smallBodies.values())+list(self.planets.values()):
            body.pos += body.velocity*self.dt

def doCollide(obj1, obj2):
    """Determine if two objects have collided."""
    # both objects have to be spheres
    return obj1.trueRadius > 0 and obj2.trueRadius > 0 and mag(obj1.pos - obj2.pos) <= (obj1.trueRadius+obj2.trueRadius)

def inelasticCollision(obj1, obj2):
    """obj 1 is larger."""
    # ghosts don't feel this force
    if obj1.ghost or obj2.ghost:
        return
    vec = (obj1.velocity - obj2.velocity)
    obj1.velocity = (obj2.mass*obj2.velocity+obj1.mass*obj1.velocity)/(obj2.mass+obj1.mass)
    obj2.visible = False

def yearsToSeconds(years):
    return years*525600*60

def vecToList(v):
    return [v.x, v.y, v.z]



######## DEFINE RUN PARAMETERS

# body initial conditions
# sun
mSun = 1.989e30
rSun = 695.508e6

# earth
mEarth = 5.972e24
aEarth = 149.6e9
eEarth = 0.0167086
rPEarth = aEarth*(1-eEarth)
rEarth = 6371e3

# apophis
mApophis = 27e9
rApophis = 2e2

# bennu

scene = canvas(title = "Asteroid Deflection", width=800, height=640, range=rPEarth*1.5)

# make SolarSystem
sun = Body('sun', mSun, rSun, rSun*9)
sun.color = color.yellow
home = SolarSystem(sun)
# make earth (body to be hit)
earth = Planet('earth', mEarth, rEarth)
earth.color = color.blue
home.addBody(earth, (0, sqrt(home.star.mu/rPEarth), 0), (rPEarth, 0, 0))
# make smallBodies
# smallBody
# def __init__(self, name, mass, trueRadius)
killer = smallBody('killer', mApophis, rApophis)
killer.color = color.red
killer.nonImpulsiveThrust = 0
ghost = smallBody('ghost', 488e9, rApophis)
ghost.color = color.white
ghost.ghost = True
# def generateInterceptor(self, body, bodyToIntercept, initPosition, tImpact)
home.generateInterceptor(killer, earth, [-earth.pos.x*3, earth.pos.x*3, 0], 1.1, 0)
# home.generateInterceptor(ghost, earth, [earth.pos.x*1.1, -earth.pos.x/3, 0], .3, 0)
savior = smallBody('savior', 1e5, 1e8)
savior.color = color.green
# home.generateInterceptor(savior, killer, vecToList(earth.pos), .8, 0)
# print(mag(home.getBody('savior').velocity-home.getBody('earth').velocity)-8000)

######## END DEFINE RUN PARAMETERS




# main loop
tvec = []
killerFromEarth = []
# a more powerful program would check every planet against every small body by making a member
# doCollide, but this is straining enough as it is
home.dt = 40
while not doCollide(home.getBody('earth'), home.getBody('killer')) and home.t < yearsToSeconds(4):
    rate(1000000)  
    home.getAccelerations()
    home.updateVelocities()
    home.updatePositions()
    killerFromEarth.append(mag(home.getBody('killer').pos - home.getBody('earth').pos) - (home.getBody('killer').trueRadius + home.getBody('earth').trueRadius))
    tvec.append(home.t/(525600*60))
    home.t += home.dt

    # years = float(home.t)/float(525600*60)
if not home.t >= yearsToSeconds(3):
    L = label(pos=vector(rPEarth*1.2, rPEarth*1.2, 0),
        text=('Stephen Hawking was right'), space=30,
        height=16,
        font='sans')
    debris = []
    body = home.getBody('earth')
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
        debris[-1].color = vector(random.random(), random.random(), random.random())

    t = 0
    dt = .001
    while t < 100:
        rate(15000)
        for deb in debris:
            deb.pos += deb.velocity*home.dt
        t += dt
plt.plot(tvec, killerFromEarth, )
plt.plot(tvec, [rEarth]*len(tvec))
plt.yscale('log')
plt.title('Asteroid Distance')
plt.xlabel('time (years)')
plt.ylabel('distance from asteroid to center of earth (m)')
plt.show()

"""Now what do we want to know? What useful data can we extract from a run of this sim?

How much did we nudge the asteroid off course at its closest point to Earth?

Can we show a 'ghost' of what would've happened had we not intervened? DONE

Graph for any given initial conditions of how much impulse we need to impart to the asteroid to avoid collision (for kinetic impactor)
create conditions, solve lambert, iterate over various impulses in the vf direction, solve tof, and get distance away from center of earth

Similar graph for non impulsive is really hard no thank you

start asteroid in some position

generate pork chop plots for departure dates

choose a good one (high vf, low energy) and launch then


"""























