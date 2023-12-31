# physics-fun
A little project to build a physics engine and use it in a small demo game

# Todo

Phase 1 (Shapes in Box);
1. [Maybe] Seperate collsion detection and handling from the sprite/initialization classes so that W_Objects have a collider object as a property rather than directly holding the collsions.
2. Create Rectangle W_object and test it's collision handling.
3. Replace Rectangle W_object with Conv_poly object

Phase 2 (Good Shapes in Box);
1. If Phase 0 - 6 or Phase 1 - 1 have not been previously completed consider doing them now, additionaly Gauss integrals for better subframe precision?
2. Implement multiple types of object materials which give different collsion properties (elastic, semi elastic, force killing different weights etc.)
3. Create Test Suite for all previously created work
4. Implement External Forces, starting with gravity, maybe also Internal Value, Coordinate and Velocity dependant forces ie Electromagnetic
5. Make one of the Shapes Controllable to provide 'player' interaction

Phase 3-K (Basic Game);
[Too far in the future and too ambitous to break down now, essentialy just create basic platformer using above tools]

Phase K+1 (3d Spheres);
[Create 3d spheres in box simulation using what I've Learned]

# Completed
At the present version, (Maybe I'll version number, in which case V.0.9.1) I have around 400 lines of code in python using the pysdl2 library. Of these lines around 150 are from early testing before I started using pysdl2 and was using trtle for graphics output {<3 Turtle graphics my love}. I am using pysdl2 for input/output handling.

The code I have so far creates a small window with boundaries placed at the edges and allows multiple balls of varying weight and radii to bounce around off eachother in fully elastic collisions.

# Notes on the Project workflow
I am working in VSCode on windows 10 using conda (miniconda3) to manage my environment, and use primarily the VSCode git interface to manage the project with occasional tweaks through the github website :smile:. After having finshed the absolute basics of the code I now want to focus on better Version Control and branching.

# Previous Phases

Phase 0 (Balls in Box); COMPLETE
1. ~~Add proper comments to code already present~~
2. ~~Reevaluate the necessity of the CVector Vector handling~~ [Removed CVector usage]
3. ~~Change refresh() actions over to using the TimeFuncs.timestep() function, as well as being tied to the system clock instead of framerate~~ [Sys Clock tied, refresh fixed in World_Objects.Ball which is the only object with proper refresh so far]
4. ~~Implement and test bsic ball on ball collision~~
5. ~~Allow free selection of Ball Size~~
6. ~~[Maybe] implement refresh() new collision rechecking to allow a Ball object to continue bouncing multiple times in one refresh~~, [Maybe] ~~implement collision based physics ticks rather than frame based ones, ie extrapolate object paths until first collision/(lazily possible collision using event horizon) and do a physics refresh exactly at that point to reduce tick frequency.~~ Second part Unfeasible at this stage
7. ~~Fix Basic Dynamic Sprite Generation for Barriers and Balls~~
8. ~~Fix physics bug (Likely in collision)~~        RESOLVED - Double counting of certain offsets was issue
9. ~~Work on gfx bugs~~
10. ~~Different weights~~
