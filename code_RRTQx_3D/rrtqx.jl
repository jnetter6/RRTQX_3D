function multirrtqx(S::Array{TS}, N::Int64, total_planning_time::Float64, slice_time::Float64,
  delta::Float64, ballConstant::Float64, changeThresh::Float64,
  searchType::String, MoveRobotFlag::Bool, saveVideoData::Bool, obstacleFile::String,
  statsArgs...) where {TS}

  T = RRTNode{Float64}

  # NOTE THIS IS HARD CODED HERE (SHOULD PROBABLY MAKE INPUT ARGUMENT)
  robotSensorRange = 10 # 20 # used for "sensing" obstacles

  if length(statsArgs) >= 2
	dataFileName = statsArgs[2]
  else
	dataFileName = "data.txt"
  end

  startTime = time_ns()
  save_elapsed_time = 0.0 # will hold save time to correct for time spent
                          # writing to files (only used for visualization video)

  ### do initialization stuff:

  # init a KDtree that will be used
  # MAKE SURE THIS USES APPROPERIATE DISTANCE FUNCTION !!!!!!!!!
  KD = []

  allDone = false

  for i = 1:N
    push!(KD, KDTree{RRTNode{Float64}}(S[i].d, KDdist))
  end
  # init queue. Note that for algorithms that do not require a queue this
  # is essentially an empty variable that is used to help with polymorphism,
  # which makes life easy, the code fast, and allows much of the same code
  # to be used for all algorithms. Note about RRTx: I decided to forgoe RRT*+
  # and focus only on #+ sice the latter will likely be more useful than the
  # former in practice

  staticObs = []
  #obstacleFile = "environments/building.txt" # Dynamic_5: ACC 20; Static_5: small envir;
  a = open(obstacleFile, "r")

  # get the number of polygons
  P = parse(Int, readline(a))

  for p = 1:P
    # get the numer of point in this polygon

    centerT = Array{Float64}(undef, 3)
    centerT[:] = str2array(readline(a))
    center = [centerT[1] centerT[2] centerT[3]]

    radius = parse(Float64, readline(a))

    obsBehavourType = parse(Int, readline(a))
    push!(staticObs, [center, radius])
  end
  close(a)

  if searchType == "RRTx"
    Q = []
    for i = 1:N
      push!(Q, rrtXQueue{RRTNode{Float64}, typeof((Float64, Float64, Float64))}())
      Q[i].Q = BinaryHeap{RRTNode{Float64}, typeof((Float64, Float64, Float64))}(keyQ, lessQ, greaterQ, markQ, unmarkQ, markedQ, setIndexQ, unsetIndexQ, getIndexQ)
      Q[i].OS = JList{RRTNode{Float64}}()
      Q[i].S = S[i]
      Q[i].changeThresh = changeThresh
    end
  else
    error("unknown search type: $(searchType)")
  end
  for i = 1:N
    S[i].sampleStack = JList{Array{Float64,3}}() # stores a stack of points that we
                                              # desire to insert in the future in
                                              # (used when an obstacle is removed
    S[i].delta = delta
  end

  robotRads = S[1].robotRadius


  # define root node in the search tree graph
  root = []
  for i = 1:N
    push!(root, RRTNode{Float64}(S[i].start))
    
    # explicit check root
    (explicitlyUnSafe, unused) = explicitNodeCheck3D(S[i], root[i])
    if explicitlyUnSafe
      error("root is not safe")
    end
  

    root[i].rrtTreeCost = 0.0
    root[i].rrtLMC = 0.0

    # insert the root into the KDtree
    kdInsert(KD[i], root[i])
  end

  goal = []

  for i = 1:N
    # define a goal node
    push!(goal, RRTNode{Float64}(S[i].goal))
    goal[i].rrtTreeCost = Inf
    goal[i].rrtLMC = Inf
    S[i].goalNode = goal[i]
    S[i].root = root[i]

    S[i].moveGoal = goal[i] # this will store a node at least as far from the root as the robot
                    # during movement it key is used to limit propogation beyond the
                    # region we care about

    S[i].moveGoal.isMoveGoal = true
  end

  moveGoals = []
  for i = 1:N
    push!(moveGoals, [])
    push!(moveGoals[i], S[i].moveGoal.position)
  end

  # paramiters that have to do with the robot path following simulation
  R = []
  for i = 1:N
    push!(R, RobotData{RRTNode{Float64}}(copy(S[i].goal), goal[i], 20000))
  end
  
  vCounter = []
  for i = 1:N
    push!(vCounter, 0) # helps with visuilizing data
    S[i].fileCtr = vCounter[i]
  end

  sliceCounter = 0 # helps with saving accurate time data

  ### end of initialization stuff


  # if saving stats about run, then allocate memory to store data
  if (length(statsArgs) >= 1 && statsArgs[1])

    savingStats = true
    estimated_number_of_data_points = 4*Int(ceil(total_planning_time/slice_time))
    checkPtr = []
    for i = 1:N
      push!(checkPtr, 1)
    end      # array location to save stats
    itOfCheck = []
    elapsedTime = []
    nodesInGraph = []
    costOfGoal = []
    for i = 1:N
      temp = Array{Int64}(undef, estimated_number_of_data_points)
      temp[1] = 0
      push!(itOfCheck, temp)

      temp = Array{Float64}(undef, estimated_number_of_data_points)
      temp[1] = 0.0
      push!(elapsedTime, temp)

      temp = Array{Int64}(undef, estimated_number_of_data_points)
      temp[1] = 1
      push!(nodesInGraph, temp)

      temp = Array{Float64}(undef, estimated_number_of_data_points)
      temp[1] = Inf
      push!(costOfGoal, temp)
    end

    #numReduces = Array{Int64}(undef, estimated_number_of_data_points)
    #numReduces[1] = 0
  else
    savingStats = false
  end

  # while planning time left, plan. (will break out when done)
  oldrrtLMC = []
  timeForGC = []
  robot_slice_start = time_ns()
  for i = 1:N
    S[i].startTimeNs = robot_slice_start
    S[i].elapsedTime = 0.0
    push!(oldrrtLMC, Inf)
    push!(timeForGC, 0)
  end

  # firstTimeQLearning = true # to check the first time to get Q_Path. default to be true

  # qPathHasCollision = false # check collision. default to be false
  for i = 1:N
    S[i].augDist = 0.0             # initialize augDist
    S[i].kino_dist = 0.0
  end

  localPoseAndKd = []
  localNormEsq = []
  localTrigCond = []
  pastVels = []
  distAndVelError = []
  allBVPs = []
  angAndTimeSince = []
  savedAngle = []
  timeSinceLastSave = []
  for i = 1:N
    push!(localPoseAndKd, Array{Tuple{Array{Float64},Float64}}(undef,1000))
    push!(localNormEsq, Array{Float64}(undef,3000))
    push!(localTrigCond, Array{Float64}(undef,3000))
    push!(pastVels, [])
    push!(distAndVelError, [])
    push!(angAndTimeSince, [])
    push!(savedAngle, 0.0)
    push!(timeSinceLastSave, 0.0)
    #push!(pastVels, Array{Tuple{Array{Float64}, Array{Float64}, Array{Float64}, Int64}}(undef, 1))
    S[i].numCoveredLocal = 0
    S[i].numLocal = 0
    S[i].numEsqTrigLocal = 0
  end
  #println(pastVels[1])
  # S.numErrTrigCoveredLocal = 0
  # S.numErrTrigLocal = 0
  ### initialize learning process


  ### end of initializing learning process

  # NormEsqvec = Array{Float64}(undef, 0)
  # TrigCondvec = Array{Float64}(undef, 0)
  lastVel = []

  for i = 1:N
    S[i].NormEsqvec = zeros(0,)
    S[i].TrigCondvec = zeros(0,)

    S[i].lastVelocity = [0.; 0.; 0.]
    push!(lastVel, [0.; 0.; 0.])
  end

  # environmentChangeFinished = false
  # while v_bot != v_goal
  hyberBallRad = []
  for i = 1:N
    push!(hyberBallRad, delta)
  end

  currPos = []
  prevPos = []
  actualPrevPos = []
  prevPosAvg = []
  for i = 1:N
    push!(currPos, R[i].robotPose)
    push!(prevPos, Array{Array{Float64}}(undef, 10))
    push!(prevPosAvg, Array{Array{Float64}}(undef, 5))
    push!(actualPrevPos, [0., 0., 0.])
    for j = 1:10
      prevPos[i][j] = currPos[i]
    end
    for j = 1:4
      prevPosAvg[i][j] = currPos[i]
    end
  end

  
  currObsPos = []
  nextObsPos = []
  nextObsPos2 = []
  nextPos = []
  nextPos2 = []
  prevPrevVel = []
  prevVel = []
  currVel = []
  gotPrevPrev = []
  gotPrev = []
  for i = 1:N
    push!(currObsPos, Array{Float64}(undef, 3, 2))
    push!(nextObsPos, Array{Float64}(undef, 3, 2))
    push!(nextObsPos2, Array{Float64}(undef, 3, 2))
    push!(prevPrevVel, Array{Float64}(undef, 3))
    push!(prevVel, Array{Float64}(undef, 3))
    push!(currVel, Array{Float64}(undef, 3))
    push!(gotPrevPrev, false)
    push!(gotPrev, false)
    push!(nextPos, R[i].robotPose)
    push!(nextPos2, R[i].robotPose)
  end

  BVPEnds = []
  maxKDs = []
  whichBVP = []
  level = []
  for i = 1:N
    tempEnds = [R[i].robotPose]
    push!(BVPEnds, tempEnds)
    push!(maxKDs, [])
    push!(maxKDs[i], 0.0)
    push!(whichBVP, 1)
    push!(level, 0)
  end

  BVPJustChanged = []
  BVPCounter = []
  NextBVPCheck = []
  for i = 1:N
    push!(BVPCounter, 0)
    push!(BVPJustChanged, false)
    push!(NextBVPCheck, true)
  end

  while true
    for i = 1:N
      if (NextBVPCheck[i] == false)
        if (BVPCounter[i] > 5)
          BVPCounter[i] = 0
          NextBVPCheck[i] = true
        else
          BVPCounter[i] += 1
        end
      end
      if (BVPJustChanged[i] == true)
        if (i == 1)
          push!(angAndTimeSince[i], [savedAngle[i], timeSinceLastSave[i], 1])
        else
          push!(angAndTimeSince[i], [savedAngle[i], timeSinceLastSave[i], 0])
        end
        timeSinceLastSave[i] = 0.0;
        push!(BVPEnds[i], currPos[i])
        BVPJustChanged[i] = false
        whichBVP[i] += 1
        push!(maxKDs[i], 0.0)
        #kdAndLV = sim_TNNLS_B_CT_Local_Max(BVPEnds[i][size(BVPEnds[i])[1] - 1][:], BVPEnds[i][size(BVPEnds[i])[1]][:], lastVel[i])
        len = dist(BVPEnds[i][size(BVPEnds[i])[1] - 1][:], BVPEnds[i][size(BVPEnds[i])[1]][:])
        #lastVel[i] = kdAndLV[2]
        println("BVP Change!")
        println(i)
        println(vCounter[i])
        #println(maxKDs[i])
        #println(length(maxKDs[i]))
        #println(whichBVP[i])
        push!(allBVPs, (len, (maxKDs[i][size(maxKDs[i])[1] - 1]), i))
        #if ((maxKDs[i][size(maxKDs[i])[1] - 1]) > (kdAndLV[1] + .1))
        #  level[i] = 1
        #  push!(allBVPs, (len, (maxKDs[i][size(maxKDs[i])[1] - 1]), 1, i))
        #else
        #  push!(allBVPs, (len, (maxKDs[i][size(maxKDs[i])[1] - 1]), 0, i))
        #end
      end
    end



    if (allDone == true)
      break
    end
    for i = 1:N
      hyberBallRad[i] = min(delta, ballConstant*((log(1+KD[i].treeSize)/(KD[i].treeSize))^(1/S[i].d)))
      itOfCheck[i][checkPtr[i]] += 1
    end
    now_time = time_ns()

    # calculate the end time of the first slice
    slice_end_time = (1+sliceCounter)*slice_time

    # see if warmup time has ended
    warmUpTimeJustEnded = false
    for i = 1:N
      if S[i].inWarmupTime && S[i].warmupTime < S[i].elapsedTime
        warmUpTimeJustEnded = true
        S[i].inWarmupTime = false
      end
    end

	### deal with larger kino_dist
  for i = 1:N
    if S[i].kino_dist > maxKDs[i][whichBVP[i]]
      maxKDs[i][whichBVP[i]] = S[i].kino_dist
    end

	  if S[i].kino_dist > S[i].augDist

      # not count in time for Minkowski sum
      before_save_time = time_ns()

      println("--------------------------------------------------------------------- New KD: $(S[i].kino_dist)")

	    S[i].augDist = S[i].kino_dist
  	  obstacleAugmentation(S[i], S[i].augDist)

	    # add latest augmented obs
  	  list_item = S[i].obstacles.front
	    while list_item != list_item.child
		  ob = list_item.data
		  if !ob.obstacleUnused
		    addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
	      end
		  list_item = list_item.child
	    end
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
      propogateDescendants(Q[i], R[i])
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")

	    if R[i].currentMoveInvalid && explicitPointCheck(S[i], R[i].robotPose)[1] # R in aug obs && R not in ori obs
		  R[i].currentMoveInvalid = false
		  println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R[i].currentMoveInvalid) ----------- robot in aug obs")
	    end
  	  if !markedOS(S[i].moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
		  verifyInQueue(Q[i], S[i].moveGoal)
	    end
	    reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])

	    save_elapsed_time += (time_ns()-before_save_time)/1000000000

	  # println("obstacle enlarged")
	  end
  end

	### end of deal with larger kino_dist


    ### add/remove newly "detected" obstacles ###
    ### beginning of remove obstacle
    if vCounter[1] > 30
    for i = 1:N
      #currObsPos[i][1] = [(currPos[i][1] - 1.0), (currPos[i][2] - 1.0)]
      #currObsPos[i][2] = [(currPos[i][1] + 1.0), (currPos[i][2] - 1.0)]
      #currObsPos[i][3] = [(currPos[i][1] + 1.0), (currPos[i][2] + 1.0)]
      #currObsPos[i][4] = [(currPos[i][1] - 1.0), (currPos[i][2] + 1.0)]
      
      #nextPos[i] = [(currPos[i][1] + 2*(currPos[i][1]-prevPos[i][5][1])) (currPos[i][2] + 2*(currPos[i][2]-prevPos[i][5][2])) (currPos[i][3] + 2*(currPos[i][3]-prevPos[i][5][3]))]
      #nextPos2[i] = [(currPos[i][1] + 4*(currPos[i][1]-prevPos[i][5][1])) (currPos[i][2] + 4*(currPos[i][2]-prevPos[i][5][2])) (currPos[i][3] + 4*(currPos[i][3]-prevPos[i][5][3]))]
      nextPos[i] = [(currPos[i][1] + 1.5*(currPos[i][1]-prevPosAvg[i][3][1])) (currPos[i][2] + 1.5*(currPos[i][2]-prevPosAvg[i][3][2])) (currPos[i][3] + 1.5*(currPos[i][3]-prevPosAvg[i][3][3]))]
      nextPos2[i] = [(currPos[i][1] + 3*(currPos[i][1]-prevPosAvg[i][3][1])) (currPos[i][2] + 3*(currPos[i][2]-prevPosAvg[i][3][2])) (currPos[i][3] + 3*(currPos[i][3]-prevPosAvg[i][3][3]))]

      #nextObsPos[i][1,:] = [(nextPos[i][1] - 1.4), (nextPos[i][2] - 1.4)]
      #nextObsPos[i][2,:] = [(nextPos[i][1] + 1.4), (nextPos[i][2] - 1.4)]
      #nextObsPos[i][3,:] = [(nextPos[i][1] + 1.4), (nextPos[i][2] + 1.4)]
      #nextObsPos[i][4,:] = [(nextPos[i][1] - 1.4), (nextPos[i][2] + 1.4)]

      #nextObsPos2[i][1,:] = [(nextPos2[i][1] - 1.4), (nextPos2[i][2] - 1.4)]
      #nextObsPos2[i][2,:] = [(nextPos2[i][1] + 1.4), (nextPos2[i][2] - 1.4)]
      #nextObsPos2[i][3,:] = [(nextPos2[i][1] + 1.4), (nextPos2[i][2] + 1.4)]
      #nextObsPos2[i][4,:] = [(nextPos2[i][1] - 1.4), (nextPos2[i][2] + 1.4)]
      currObs = SphereObstacle(currPos[i], (1.5))
      nextObs = SphereObstacle(nextPos[i], (3.0))
      nextObs2 = SphereObstacle(nextPos2[i], (3.0))
    
      currObs.startTime = S[1].elapsedTime
      currObs.lifeSpan = slice_time*3
      currObs.obstacleUnused = false

      nextObs.startTime = S[1].elapsedTime
      nextObs.lifeSpan = slice_time*3
      nextObs.obstacleUnused = false

      nextObs2.startTime = S[1].elapsedTime
      nextObs2.lifeSpan = slice_time*3
      nextObs2.obstacleUnused = false

      #if (i != 2)
      #  if (Wdist(R[2].robotPose, currPos[i]) < robotSensorRange)
      #    addObsToCSpace(S[2], currObs)
      #    if (Wdist(R[2].robotPose, nextPos[i]) > 1.5)
      #      addObsToCSpace(S[2], nextObs)
      #    end
      #    if (Wdist(R[2].robotPose, nextPos2[i]) > 1.5)
      #      addObsToCSpace(S[2], nextObs2)
      #    end
      #  end
      #end
      if (i != 1)
        if (Wdist(R[1].robotPose, currPos[i]) < robotSensorRange)
          #if (level[i] == 0)
            addObsToCSpace(S[1], currObs)
          #end
          if (Wdist(R[1].robotPose, nextPos[i]) > 1.5)
            #if (level[i] == 0)
              addObsToCSpace(S[1], nextObs)
            #end
          end
          if (Wdist(R[1].robotPose, nextPos2[i]) > 1.5)
            #if (level[i] == 0)
              addObsToCSpace(S[1], nextObs2)
            #end
          end
        end
      end
      if (i != 2)
        if (Wdist(R[2].robotPose, currPos[i]) < robotSensorRange)
          #if (level[i] == 0)
            addObsToCSpace(S[2], currObs)
          #end
          if (Wdist(R[2].robotPose, nextPos[i]) > 1.5)
            #if (level[i] == 0)
              addObsToCSpace(S[2], nextObs)
            #end
          end
          if (Wdist(R[2].robotPose, nextPos2[i]) > 1.5)
            #if (level[i] == 0)
              addObsToCSpace(S[2], nextObs2)
            #end
          end
        end
      end
    end
    end
    # remove obstacles at the required time
    for i = 1:N
    S[i].elapsedTime = (time_ns() - S[i].startTimeNs)/1000000000 - save_elapsed_time
    list_item = S[i].obstacles.front
    removedObstacle = false
    while list_item != list_item.child
      ob = list_item.data

      if !ob.senseableObstacle && !ob.obstacleUnused && (ob.startTime + ob.lifeSpan <= S[i].elapsedTime)
        # time to remove obstacle
        removeObstacle(S[i], KD[i], Q[i], ob, root[i], hyberBallRad[i], S[i].elapsedTime, S[i].moveGoal)
        removedObstacle = true
      elseif ob.senseableObstacle && ob.obstacleUnusedAfterSense && Wdist(R[i].robotPose, ob.position) < robotSensorRange + ob.radius
        # place to remove obstacle

        # because the space that used to be in this obstacle was never sampled
        # there will be a hole in the graph where it used to be. The following
        # attempts to mitigate this problem by requiring that the next few samples
        # come from the space that used to be inside the obstacle
        randomSampleObs(S[i], KD[i], ob) # stores samples in the sample stack
        removeObstacle(S[i], KD[i], Q[i], ob, root[i], hyberBallRad[i], S[i].elapsedTime, S[i].moveGoal)
        ob.senseableObstacle = false
        ob.startTime = Inf
        removedObstacle = true
      elseif S[i].spaceHasTime && ob.nextDirectionChangeTime > R[i].robotPose[4] && ob.lastDirectionChangeTime != R[i].robotPose[4]
        # a moving obstacle with unknown path is changing direction, so remove
        # its old anticipated trajectory

        removeObstacle(S[i], KD[i], Q[i], ob, root[i], hyberBallRad[i], S[i].elapsedTime, S[i].moveGoal)
        ob.obstacleUnused = false # obstacle is still used
        removedObstacle = true
      end

      list_item = list_item.child
    end
  

	  # if S.elapsedTime >= 13.0 && !environmentChangeFinished
	  # 	ob = S.obstacles.front.data
	  # 	randomSampleObs(S, KD, ob)
	  # 	removeObstacle(S, KD, Q, ob, root, hyberBallRad, S.elapsedTime, S.moveGoal)
	  # 	removedObstacle = true
	  # 	environmentChangeFinished = true
	  # end

	  if removedObstacle
      #println("----------------------------------------------------------------------------- Removed obstacle")
      reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
    end
  end
    ### end of remove obstacle

    ### beginning of add obstacle
    # add obstacles at the required time
  for i = 1:N
    list_item = S[i].obstacles.front
    addedObstacle = false
    while list_item != list_item.child
      ob = list_item.data
      if !ob.senseableObstacle && ob.obstacleUnused && (ob.startTime <= S[i].elapsedTime <= ob.startTime + ob.lifeSpan)
        # time to add
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        addedObstacle = true
      elseif ob.senseableObstacle && !ob.obstacleUnusedAfterSense && Wdist(R[i].robotPose, ob.position) < robotSensorRange + ob.radius
        # place to add obstacle
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        ob.senseableObstacle = false
        addedObstacle = true
      elseif S[i].spaceHasTime && ob.nextDirectionChangeTime > R[i].robotPose[3] && ob.lastDirectionChangeTime != R[i].robotPose[3]
        # time that a moving obstacle with unknown path changes direction
        ob.obstacleUnused = false
        changeObstacleDirection(S[i], ob, R[i].robotPose[3])
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        ob.lastDirectionChangeTime = copy(R[i].robotPose[3])
        #println("$(ob.nextDirectionChangeTime)  $(S.moveGoal.position[3]) ")
        addedObstacle = true
      elseif warmUpTimeJustEnded && !ob.obstacleUnused
        # warm up time is over, so we need to treat all active obstacles
        # as if they have just been added
        addNewObstacle(S[i], KD[i], Q[i], ob, root[i], vCounter[i], R[i])
        addedObstacle = true
      end

      list_item = list_item.child
    end
    if addedObstacle
      # propogate inf cost to all nodes beyond the obstacle and in its
      # basin of attraction
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
	  propogateDescendants(Q[i], R[i])
	  # println("-------------------------------------------------------------------------------R.currentMoveInvalid = $(R.currentMoveInvalid)")
      if !markedOS(S[i].moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
        verifyInQueue(Q[i], S[i].moveGoal)
      end
      #println("--------------------------------------------------------------------------------- Added obstacle")
      reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
    end
  end
    ### end of add obstacle

# 	if (S.kino_dist > augDist) || removedObstacle || addedObstacle
#
# 	  obstacleAugmented = false
#
# 	  if S.kino_dist > augDist
# 		println("kd increased")
# 		augDist = S.kino_dist
# 	  end
#
# 	  # pop augmented obs if needed
# 	  listEmpty(S.augObs)
# 	  # now S.obstacles only contains original obs, S.augObs is empty
# ##################################################################################################
# 	  # push augmented obs into S.obstacles and S.augObs
# 	  obstacleAugmentation(S, augDist)
#
# 	  list_item = S.augObs.front
# 	  while list_item != list_item.child
# 		  ob = list_item.data
# 		  addNewObstacle(S, KD, Q, ob, root, vCounter[i], R)
# 		  list_item = list_item.child
# 	  end
# 	  propogateDescendants(Q, R)
# 	  if !markedOS(S.moveGoal) # I'm pretty sure this is always true, since OS is emopty here -- M.O.
# 		  verifyInQueue(Q, S.moveGoal)
# 	  end
# 	  reduceInconsistency(Q, S.moveGoal, robotRads, root, hyberBallRad)
# 	  obstacleAugmented = true
# 	  println("added augmented obstacles")
#
#     end
    ### done with add/remove newly "detected" obstacles ###

    allDone = true
    for i = 1:N
      if R[i].robotPose == root[i].position
        break
      else
        allDone = false
      end
    end

    for i = 1:N
    # if this robot has used all of its allotted planning time of this slice
    S[i].elapsedTime = (time_ns() - S[i].startTimeNs)/1000000000 - save_elapsed_time
    if (S[i].elapsedTime >= slice_end_time)

      # calculate the end time of the next slice
      slice_end_time = (1+sliceCounter)*slice_time

      robot_slice_start = now_time
      if i == 1
        sliceCounter += 1
      end

      truncElapsedTime = floor(S[i].elapsedTime * 1000)/1000
      if i == 1
        println("slice $(sliceCounter) --- $(truncElapsedTime) -------- $(S[i].moveGoal.rrtTreeCost) $(S[i].moveGoal.rrtLMC) ----")
      end

      for j = 10:-1:2
        prevPos[i][j] = prevPos[i][j-1]
      end
      currPos[i] = R[i].robotPose
      prevPos[i][1] = currPos[i]
      
      if (i > 1)
        if (actualPrevPos[i] != currPos[i])
          playerAgentDist = sqrt((prevPos[i][1][1] - prevPos[1][1][1])^2 + (prevPos[i][1][2] - prevPos[1][1][2])^2 + (prevPos[i][1][3] - prevPos[1][1][3])^2)
          #println(playerAgentDist)
          prevPos[i][1] = [(prevPos[i][1][1] + (0.01*rand(Float64) - .005)*playerAgentDist), (prevPos[i][1][2] + (0.01*rand(Float64) - .005)*playerAgentDist), (prevPos[i][1][3] + (0.01*rand(Float64) - .005)*playerAgentDist)]
        end
      end

      for j = 5:-1:1
        prevPosAvg[i][j] = [(((prevPos[i][(2*j)][1]) + (prevPos[i][(2*j-1)][1]))/2),(((prevPos[i][(2*j)][2]) + (prevPos[i][(2*j-1)][2]))/2), (((prevPos[i][(2*j)][3]) + (prevPos[i][(2*j-1)][3]))/2)]
      end

      actualPrevPos[i] = currPos[i]

      if (vCounter[i] > 35)
        #playerAgentDist = sqrt((prevPos[i][1][1] - prevPos[1][1][1])^2 + (prevPos[i][1][2] - prevPos[1][1][2])^2 + (prevPos[i][1][3] - prevPos[1][1][3])^2)
        #noisyPastVec1 = [(prevPos[i][4][1] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][4][2] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][4][3] + (0.02*rand(Float64) - .01)*playerAgentDist)]
        #noisyPastVec2 = [(prevPos[i][5][1] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][5][2] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][5][3] + (0.02*rand(Float64) - .01)*playerAgentDist)]
        #noisyCurrVec1 = [(prevPos[i][1][1] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][1][2] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][1][3] + (0.02*rand(Float64) - .01)*playerAgentDist)]
        #noisyCurrVec2 = [(prevPos[i][2][1] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][2][2] + (0.02*rand(Float64) - .01)*playerAgentDist), (prevPos[i][2][3] + (0.02*rand(Float64) - .01)*playerAgentDist)]
        #pastVec = [(prevPos[i][4][1] - prevPos[i][5][1]), (prevPos[i][4][2] - prevPos[i][5][2]), (prevPos[i][4][3] - prevPos[i][5][3])]
        #currVec = [(prevPos[i][1][1] - prevPos[i][2][1]), (prevPos[i][1][2] - prevPos[i][2][2]), (prevPos[i][1][3] - prevPos[i][2][3])]
        pastVec = [(prevPosAvg[i][4][1] - prevPosAvg[i][5][1]), (prevPosAvg[i][4][2] - prevPosAvg[i][5][2]), (prevPosAvg[i][4][3] - prevPosAvg[i][5][3])]
        currVec = [(prevPosAvg[i][1][1] - prevPosAvg[i][2][1]), (prevPosAvg[i][1][2] - prevPosAvg[i][2][2]), (prevPosAvg[i][1][3] - prevPosAvg[i][2][3])]
        #pastVec = [(noisyPastVec1[1] - noisyPastVec2[1]), (noisyPastVec1[2] - noisyPastVec2[2]), (noisyPastVec1[3] - noisyPastVec2[3])]
        #currVec = [(noisyCurrVec1[1] - noisyCurrVec2[1]), (noisyCurrVec1[2] - noisyCurrVec2[2]), (noisyCurrVec1[3] - noisyCurrVec2[3])]
        angle = acos(((pastVec[1]*currVec[1]) + (pastVec[2]*currVec[2]) + (pastVec[3]*currVec[3]))/(sqrt(pastVec[1]^2 + pastVec[2]^2 + pastVec[3]^2)*sqrt(currVec[1]^2 + currVec[2]^2 + currVec[3]^2)))
        angle = angle*(360/(2*pi))
        if ((abs(angle) > 18) && (NextBVPCheck[i] == true))
          BVPJustChanged[i] = true
          NextBVPCheck[i] = false
          savedAngle[i] = angle
        else
          timeSinceLastSave[i] = timeSinceLastSave[i] + 1.0
        end
        if (mod(vCounter[i], 5) == 0)
          if (gotPrevPrev[i] == false)
            gotPrevPrev[i] = true
            prevPrevVel[i] = currVec
          elseif (gotPrev[i] == false)
            gotPrev[i] = true
            prevVel[i] = currVec
          else
            currVel[i] = currVec
            #estPos = [(prevPos[i][5][1] + (prevVel[i][1] * 5 * slice_time)), (prevPos[i][5][2] + (prevVel[i][2] * 5 * slice_time)), (prevPos[i][5][3] + (prevVel[i][3] * 5 * slice_time))]
            estPos = [(prevPosAvg[i][3][1] + (prevVel[i][1] * 5 * slice_time)), (prevPosAvg[i][3][2] + (prevVel[i][2] * 5 * slice_time)), (prevPosAvg[i][3][3] + (prevVel[i][3] * 5 * slice_time))]
            posError = sqrt((estPos[1] - currPos[i][1])^2 + (estPos[2] - currPos[i][2])^2 + (estPos[3] - currPos[i][3])^2)
            if (i == 1)
              push!(pastVels[i], (prevPrevVel[i], prevVel[i], currVel[i], 1))
              push!(distAndVelError[i], (posError, norm(prevVel[i]), 1))
            else
              push!(pastVels[i], (prevPrevVel[i], prevVel[i], currVel[i], 0))
              push!(distAndVelError[i], (posError, norm(prevVel[i]), 0))
            end
            prevPrevVel[i] = prevVel[i]
            prevVel[i] = currVel[i]
            #println(pastVels[i])
            #println(typeof(pastVels[i]))
            #println(typeof(pastVels[i][1]))
            #println(pastVels[i][1])
          end
        end
      end

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # record data
        elapsedTime[i][checkPtr[i]] = S[i].elapsedTime
      end

      ## move robot if the robot is allowed to move, otherwise planning is finished
      # so break out of the control loop
      if elapsedTime[i][checkPtr[i]] > total_planning_time + slice_time
        if MoveRobotFlag
          moveRobot_Q(S[i], Q[i], KD[i], slice_time, root[i], hyberBallRad[i], R[i], localPoseAndKd[i], localNormEsq[i], localTrigCond[i],save_elapsed_time)#, NormEsqvec, TrigCondvec) # 2 steps, update S.kino_dist
        else
          println("done (not moving robot)")
          break
        end
      end

      if (last(moveGoals[i]) != S[i].moveGoal.position)
        push!(moveGoals[i], S[i].moveGoal.position)
      end

      if searchType == "RRT#" || searchType == "RRTx"
        reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
        if (S[i].moveGoal.rrtLMC != oldrrtLMC[i])
          oldrrtLMC[i] = (S[i].moveGoal.rrtLMC)
        end
	    end

      ## visualize graph #############
      if saveVideoData
        before_save_time = time_ns()

		# for visualization, S.obstacles only contains original obs
		# if obstacleAugmented
		# 	saveObstacleLocations(S.augObs, "temp/augObs_$(vCounter[i]).txt")
		# 	for i = 1:S.augObs.length
		# 		listPop(S.obstacles)
		# 	end
		# end

        saveRRTTree(KD[i], "temp/edges_$(i)_$(vCounter[i]).txt")
        saveRRTNodes(KD[i], "temp/nodes_$(i)_$(vCounter[i]).txt")
        #saveRRTNodesCollision(KD, "temp/Cnodes_$(vCounter[i]).txt")
        saveRRTPath_Q(S[i], S[i].moveGoal, root[i], R[i], "temp/path_$(i)_$(vCounter[i]).txt")
        saveObstacleLocations(S[i].obstacles, "temp/obstacles_$(i)_$(vCounter[i]).txt")
		    saveOriginalObstacleLocations_Q(S[i].obstacles, "temp/originalObs_$(i)_$(vCounter[i]).txt")
        saveData(R[i].robotMovePath[1:R[i].numRobotMovePoints,:], "temp/robotMovePath_$(i)_$(vCounter[i]).txt")
		    saveKds_Q(S[i], "temp/kd_$(i)_$(vCounter[i]).txt")

        
        vCounter[i] += 1
        S[i].fileCtr = vCounter[i]
        vCounter[i] = vCounter[1]

        save_elapsed_time += (time_ns()-before_save_time)/1000000000
      end
      ## end of visualize graph ######

      # check if the robot has reached its movement goal

      # if saving stats
      if length(statsArgs) >= 1 && statsArgs[1]
        # update statistics about run, assuming that we are saving them

        if checkPtr[i] < length(costOfGoal[i])
          checkPtr[i] += 1
          itOfCheck[i][checkPtr[i]] = itOfCheck[i][(checkPtr[i]-1)] + 1

          nodesInGraph[i][checkPtr[i]] = KD[i].treeSize
          costOfGoal[i][checkPtr[i]] = min(goal[i].rrtTreeCost, goal[i].rrtLMC)
          #costOfGoal[checkPtr[i]] = extractPathLength(goal , root)
          #numReduces[checkPtr[i]] = Q.numReduces
        else
          #println("WARNING: out of space to save stats")
        end
      end
    end
  end

    #### END of obstacle and robot pose update
    #### START of normal graph search stuff

    for i = 1:N
    # pick a random node
    newNode = S[i].randNode(S[i])

    if newNode.kdInTree # happens when we explicitly sample the goal every so often
      # nodes will be updated automatically as information gets propogated to it
      continue
    end


    # find closest old node to the new node
    (closestNode, closestDist) = kdFindNearest(KD[i], newNode.position)

    # saturate
    #if closestDist > delta && newNode != S.goalNode
    #  newNode.position = closestNode.position  + (newNode.position - closestNode.position)*delta/closestDist
    #end

    if closestDist > delta && newNode != S[i].goalNode
      saturate(newNode.position, closestNode.position, delta)
    end



    # check for collisions vs static obstacles
    (explicitlyUnSafe, retCert) = explicitNodeCheck(S[i], newNode)

    if explicitlyUnSafe
      continue
    end

    #!!! Need to look into this
    GC.enable(false)

    # extend
    extend(S[i], KD[i], Q[i], newNode, closestNode, delta, hyberBallRad[i], S[i].moveGoal)



    # make graph consistant (RRT# and RRTx)
    if searchType == "RRT#" || searchType == "RRTx"
      reduceInconsistency(Q[i], S[i].moveGoal, robotRads, root[i], hyberBallRad[i])
      if(S[i].moveGoal.rrtLMC != oldrrtLMC[i])
        #printPathLengths(S.moveGoal)
        oldrrtLMC[i] = S[i].moveGoal.rrtLMC
      end
    end

    GC.enable(true)
    end
  end

  ## end of while(true)
  for i = 1:N
  elapsedTime[i][checkPtr[i]] = (time_ns()-startTime)/1000000000

  if (length(statsArgs) >= 1 && statsArgs[1])
    if (!goal[i].rrtParentUsed)
      print("goal has no parent")
    end

    stats = hcat(elapsedTime[i], itOfCheck[i], nodesInGraph[i], costOfGoal[i])

    #saveData(stats[1:checkPtr[i],:], dataFileName[i])

    #reduceData = [numReduces', nodesInGraph']'
    #saveData(reduceData[1:checkPtr,:], "temp/reduceStats.txt")
  end
  moveLength = 0
  moveLength = sum(sqrt, sum((R[i].robotMovePath[1:R[i].numRobotMovePoints-1, :] - R[i].robotMovePath[2:R[i].numRobotMovePoints, :]).^2, dims=2))

  println("distance traveled by robot: $(moveLength[1])")
  println("KD_max: $(S[i].augDist)")
  end
  #println(distAndVelError[1])
  for i = 1:N
    #saveVels(pastVels[i], "temp/AApastVels_$(i).txt")
    #saveErrors(distAndVelError[i], "temp/AAErrors_$(i).txt")
    saveErrors(angAndTimeSince[i], "temp/angles_$(i).txt")
    distances = findClosestObs(BVPEnds[i], staticObs)
    saveBVPEnds(BVPEnds[i], "temp/BVPEnds_$(i).txt")
    saveBVPDists(distances, "temp/BVPDistFromStaticObs_$(i).txt")
    saveMoveGoals(moveGoals[i], "temp/MoveGoals_$(i).txt")
  end
  saveBVPs(allBVPs, "temp/BVPs.txt")
  #println(BVPEnds[1])
  #println(maxKDs)
  #println(level)
  #println(pastVels[1])
  return (S[1].NormEsqvec, S[1].TrigCondvec)
  # saveData(tr, "temp/Trig.txt")
  # saveData(er, "temp/Esq.txt")
end
