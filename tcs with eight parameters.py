from SALib.sample.saltelli import sample
from SALib.analyze import sobol
import numpy as np
import math
import tcsode as simple_model
import timeit
import random
from multiprocessing import Process, Queue, current_process, freeze_support

PROC_NUM = 4
samples_to_take = 512
name = 'tcs'


def main(name):
    name1 = name.lower()
    parameters, bounds = setups(name1)
    problem = {'num_vars': len(parameters),
               'names': parameters,
               'bounds': bounds}

    params = sample(problem, samples_to_take, calc_second_order=True)

    print(len(params))
    if False:
        params = log_distribute(params)

    done = Queue()
    paramsdiv = []
    p = []
    remaining = 0
    for i in range(0, PROC_NUM):
        x = params.shape[0] / PROC_NUM
        remaining += x - int(x)
        x = int(x)
        y = params[x * i:x * (i + 1), :]

        if i == PROC_NUM - 1 and remaining != 0:
            y = np.append(y, params[-int(remaining):, :], axis=0)
        paramsdiv.append(y)
        p.append(Process(target=evaluate, args=(paramsdiv[i], done, i)))
        p[i].start()
    print(np.asarray(paramsdiv).shape)

    Ydiv = np.zeros(PROC_NUM).tolist()
    count = 0
    while True:
        temp, no = done.get()
        Ydiv[no] = temp.tolist()
        count += 1
        if count == PROC_NUM:
            break
    Y = []
    for i in range(0, PROC_NUM):
        Y += Ydiv[i]

    Y = np.asarray(Y)
    print(Y.shape)

    Si = sobol.analyze(problem, Y, print_to_console=True, calc_second_order=True)
    print(Si)
    Si.plot()
    write_file(Si, parameters, name)


def log_distribute(params):
    new = []
    for set in params:
        temp = []
        for i, x in enumerate(set):
            temp.append(math.log(x, 10))
        new.append(temp)

    mean = []
    std = []
    for i in range(len(params[0])):
        mean.append(np.mean(new[:][i]))
        std.append(np.std(new[:][i]))

    for i, x in enumerate(new):
        for j, y in enumerate(x):
            Z = random.gauss(0, 1)
            new[i][j] = 10 ** (mean[j] + std[j] * Z)
    return np.asarray(new)


def evaluate(values, done, no):
    Y = np.zeros([values.shape[0]])
    print("Start evaluation method with " + str(len(values)) + " parameter sets")
    start = timeit.default_timer()
    # guess = [100,100,100,100,100,100,100,100,100,100,100,100,100,100,100,100]
    for i, X in enumerate(values):
        Y[i] = simple_model.main(X[0], X[1], X[2], X[3], X[4], X[5], X[6], X[7])
        # guess = Y[i]
        if i % 100 == 0 and i != 0:
            time_100 = timeit.default_timer() - start
            remaining_time = (int((len(values) - i) / 100) * time_100) / 60
            print('progress: ' + str((i / float(len(values))) * 100) + ' %')
            print('time for 100 calculations: ' + str(time_100) + ' sec')
            print('remaining time: ' + str(remaining_time) + ' min\n')
            start = timeit.default_timer()
    done.put([Y, no])


def write_file(Si, parameters, name):
    filename = 'tcs_512_8_3.csv'
    columns = ("Param1,Param2,S1,S1con,S2,S2con,ST,STcon")
    file = open(filename, 'w')
    file.write(columns + '\n')
    S1 = Si['S1']
    S1con = Si['S1_conf']
    ST = Si['ST']
    STcon = Si['ST_conf']
    S2 = Si['S2'].tolist()
    S2con = Si['S2_conf'].tolist()
    for i, x in enumerate(S1):
        file.write(parameters[i] + ', - ,' + str(S1[i]) + ',' + str(S1con[i]) + ', - , - ,' + str(ST[i]) + ',' + str(
            STcon[i]) + '\n')
    for i in range(0, len(parameters)):
        for j in range(0, len(parameters)):
            if math.isnan(S2[i][j]) == False:
                file.write(parameters[i] + ',' + parameters[j] + ', - , - ,' + str(S2[i][j]) + ',' + str(
                    S2con[i][j]) + ', - , - ,' + '\n')
    file.close()


def setups(conf):
    parameters = ['kd2', 'kb2', 'kph', 'kd3', 'kb3', 'kbLZ', 'kdLZ', 'kbSH3']
    b = [[1, 10],
         [0.1, 1],
         [0.1, 1],
         [0.1, 1],
         [0.1, 0.3],
         [0.1, 1],
         [0.1, 1],
         [0.1, 0.3]]
    return parameters, b


if __name__ == '__main__':
    main(name)