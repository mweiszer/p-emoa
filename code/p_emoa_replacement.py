#this function implements replacement procedure for Inspyred 1.0.1 package for Python: https://inspyred.github.io See replacers.py in the package for more info.

def p_emoa_replacement(random, population, parents, offspring, args):
    """Replaces population using the biased non-dominated sorting technique.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
       weight_vector -- the list of weights 
    
    """
    weight_vector = args.setdefault('weight_vector', [0.469, 0.71])
    upper_cost_bounds = [1.2*u for u in weight_vector]
    lower_cost_bounds = [0.8*u for u in weight_vector]
    weight_vector_bounds = []
    zipped = zip(upper_cost_bounds,lower_cost_bounds)
    for cost_index1 in range(len(upper_cost_bounds)):
        for cost_index2 in range(len(lower_cost_bounds)):
            weight_vector_bounds.append([zipped[0][cost_index1],zipped[1][cost_index2]])
    
    survivors = []
    combined = list(population)
    combined.extend(offspring)
    
    # Perform the non-dominated sorting to determine the fronts.
    fronts = []
    pop = set(range(len(combined)))
    while len(pop) > 0:
        front = []
        for p in pop:
            dominated = False
            for q in pop:
                if combined[p] < combined[q]:
                    dominated = True
                    break
            if not dominated:
                front.append(p)
        fronts.append([dict(individual=combined[f], index=f) for f in front])
        pop = pop - set(front)
    
    # Go through each front and add all the elements until doing so
    # would put you above the population limit. At that point, fall
    # back to the crowding distance to determine who to put into the
    # next population. 
    for front_i, front in enumerate(fronts):
        if True:
            # Determine the crowding distance.
            distance = [0 for _ in range(len(combined))]
            individuals = list(front)
            num_individuals = len(individuals)
            num_objectives = len(individuals[0]['individual'].fitness)
            for obj in range(num_objectives):
                individuals.sort(key=lambda x: x['individual'].fitness[obj])
                distance[individuals[0]['index']] = float('inf')
                distance[individuals[-1]['index']] = float('inf')
                for i in range(1, num_individuals-1):
		    if (float(individuals[-1]['individual'].fitness[obj]-individuals[0]['individual'].fitness[obj]))<> 0:
			distance[individuals[i]['index']] = (distance[individuals[i]['index']] + 
                                                         (individuals[i+1]['individual'].fitness[obj] - 
                                                          individuals[i-1]['individual'].fitness[obj])/(float(individuals[-1]['individual'].fitness[obj]-individuals[0]['individual'].fitness[obj])))
                    else:
			distance[individuals[i]['index']] = (distance[individuals[i]['index']] + 
                                                         (individuals[i+1]['individual'].fitness[obj] - 
                                                          individuals[i-1]['individual'].fitness[obj])/(float(individuals[-1]['individual'].fitness[obj])))
            individuals.sort(key=lambda x: x['individual'].fitness[0]*weight_vector[0]+x['individual'].fitness[0]*weight_vector[1])
            cost = [sum([individuals[i]['individual'].fitness[obj]*weight_vector[obj] for obj in range(0,len(weight_vector))]) for i in range(num_individuals)]
            
            #this is currently implemented only for weights for 2 objectives
            cost1 = [sum([individuals[i]['individual'].fitness[obj]*weight_vector_bounds[0][obj] for obj in range(0,len(weight_vector_bounds[0]))]) for i in range(num_individuals)]
            cost2 = [sum([individuals[i]['individual'].fitness[obj]*weight_vector_bounds[1][obj] for obj in range(0,len(weight_vector_bounds[1]))]) for i in range(num_individuals)]
            cost3 = [sum([individuals[i]['individual'].fitness[obj]*weight_vector_bounds[2][obj] for obj in range(0,len(weight_vector_bounds[2]))]) for i in range(num_individuals)]
            cost4 = [sum([individuals[i]['individual'].fitness[obj]*weight_vector_bounds[3][obj] for obj in range(0,len(weight_vector_bounds[3]))]) for i in range(num_individuals)]
            center = cost.index(min(cost))
            bound1 = cost1.index(min(cost1))
            sol_bound1 = individuals[cost1.index(min(cost1))]
            bound2 = cost2.index(min(cost2))
            sol_bound2 = individuals[cost2.index(min(cost2))]
            bound3 = cost3.index(min(cost3))
            sol_bound3 = individuals[cost3.index(min(cost3))]
            bound4 = cost4.index(min(cost4))
            sol_bound4 = individuals[cost4.index(min(cost4))]
            f1_bound = [min([sol_bound1['individual'].fitness[0],sol_bound2['individual'].fitness[0],sol_bound3['individual'].fitness[0],sol_bound4['individual'].fitness[0]]), max([sol_bound1['individual'].fitness[0],sol_bound2['individual'].fitness[0],sol_bound3['individual'].fitness[0],sol_bound4['individual'].fitness[0]])]
            f2_bound = [min([sol_bound1['individual'].fitness[1],sol_bound2['individual'].fitness[1],sol_bound3['individual'].fitness[1],sol_bound4['individual'].fitness[1]]), max([sol_bound1['individual'].fitness[1],sol_bound2['individual'].fitness[1],sol_bound3['individual'].fitness[1],sol_bound4['individual'].fitness[1]])]
            bounds = [center]
            bounds.sort()
            
            #calculate the modified crowding distance based on weights

            for i in range(0, num_individuals):
		    if i == bounds[0]: #center point
			
			distance[individuals[i]['index']] = float('inf')
		    elif individuals[i]['individual'].fitness[0] >= f1_bound[0] and individuals[i]['individual'].fitness[0] <= f1_bound[1] and individuals[i]['individual'].fitness[1] >= f2_bound[0] and individuals[i]['individual'].fitness[1] <= f2_bound[1]:
			distance[individuals[i]['index']] = 1000000+distance[individuals[i]['index']]
		    else:
			distance[individuals[i]['index']] = 1/cost[i] #distance[individuals[i]['index']] #                
            crowd = [dict(dist=distance[f['index']], index=f['index']) for f in front]
            crowd.sort(key=lambda x: x['dist'], reverse=True)

            last_rank = [combined[c['index']] for c in crowd]
            r = 0
            num_added = 0
            num_left_to_add = len(population) - len(survivors)
            while r < len(last_rank) and num_added < num_left_to_add:
                if last_rank[r] not in survivors:
                    survivors.append(last_rank[r])
                    num_added += 1
                r += 1
            # If we've filled out our survivor list, then stop.

            if len(survivors) == len(population):
                break

    return survivors