import sys
from deap import algorithms
from deap import tools
import csv

class optimiser():

    interrupt = False

    def finish_running(self, signal, frame):
        self.interrupt = True
        print('Finishing after this generation!')

    def eaSimple(self, population, toolbox, cxpb, mutpb, ngen, stats=None,
                 halloffame=None, hoffile=None, verbose=__debug__):

        evaluationID = 0

        logbook = tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        ids = range(evaluationID, evaluationID + len(invalid_ind))
        evaluationID += len(invalid_ind)
        for ind, fit in zip(invalid_ind, fitnesses):
            ind.fitness.values = fit

        if halloffame is not None:
            halloffame.update(population)
            with open(hoffile,'wb') as out:
                csv_out=csv.writer(out)
                for row in halloffame:
                    row.append(0)
                    csv_out.writerow(row)

        record = stats.compile(population) if stats else {}
        logbook.record(gen=0, nevals=len(invalid_ind), **record)
        if verbose:
            print logbook.stream

        # Begin the generational process
        for gen in range(1, ngen + 1):
            if self.interrupt:
                self.interrupt = False
                break
            # Select the next generation individuals
            offspring = toolbox.select(population, len(population))

            # Vary the pool of individuals
            offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            ids = range(evaluationID, evaluationID + len(invalid_ind))
            evaluationID += len(invalid_ind)
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses, ids):
                ind.fitness.values = fit

            # Update the hall of fame with the generated individuals
            if halloffame is not None:
                halloffame.update(offspring)
                with open(hoffile,'a') as out:
                    csv_out=csv.writer(out)
                    for row in halloffame:
                        row.append(gen)
                        csv_out.writerow(row)

            # Replace the current population by the offspring
            population[:] = offspring

            # Append the current generation statistics to the logbook
            record = stats.compile(population) if stats else {}
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)
            if verbose:
                print logbook.stream

        return population, logbook

    def eaMuPlusLambda(self, population, toolbox, mu, lambda_, cxpb, mutpb, ngen, stats=None, halloffame=None, hoffile=None, verbose=__debug__):

        evaluationID = 0

        logbook = tools.Logbook()
        logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

        # Evaluate the individuals with an invalid fitness
        invalid_ind = [ind for ind in population if not ind.fitness.valid]
        ids = range(evaluationID, evaluationID + len(invalid_ind))
        for ind, id in zip(invalid_ind, ids):
            ind.id = id
        evaluationID += len(invalid_ind)
        fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
        for ind, fit, id in zip(invalid_ind, fitnesses, ids):
            ind.fitness.values = fit

        if halloffame is not None:
            halloffame.update(population)
            with open(hoffile,'wb') as out:
                csv_out=csv.writer(out)
                for row in halloffame:
                    row.append(0)
                    row.append(row.id)
                    csv_out.writerow(row)

        record = stats.compile(population) if stats is not None else {}
        logbook.record(gen=0, nevals=len(invalid_ind), **record)
        if verbose:
            print logbook.stream

        # Begin the generational process
        for gen in range(1, ngen + 1):
            # Vary the population
            offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)

            # Evaluate the individuals with an invalid fitness
            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            ids = range(evaluationID, evaluationID + len(invalid_ind))
            for ind, id in zip(invalid_ind, ids):
                ind.id = id
            evaluationID += len(invalid_ind)
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit, id in zip(invalid_ind, fitnesses, ids):
                ind.fitness.values = fit

            # Update the hall of fame with the generated individuals
            if halloffame is not None:
                halloffame.update(offspring)
                with open(hoffile,'a') as out:
                    csv_out=csv.writer(out)
                    for row in halloffame:
                        row.append(gen)
                        row.append(row.id)
                        csv_out.writerow(row)

            # Select the next generation population
            population[:] = toolbox.select(population + offspring, mu)

            # Update the statistics with the new population
            record = stats.compile(population) if stats is not None else {}
            logbook.record(gen=gen, nevals=len(invalid_ind), **record)
            if verbose:
                print logbook.stream

        return population, logbook
