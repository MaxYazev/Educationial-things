# -*- coding: utf-8 -*-
"""
Расчет констуркции моста, состоящего из треугольных ферм, методом конечных элементов
"""
import matplotlib.pyplot as plt
import numpy as np
import math

class Point():
    """Контейнер для точки (узлы)"""
    
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.force_x = 0
        self.force_y = 0
        self.u = 0 # Смещение по x
        self.v = 0 # Смещение по y
        self.is_fixed = False
    
    def draw(self):
        plt.plot(self.x, 
                 self.y, 
                 marker = 'o', 
                 markeredgecolor = 'k',
                 markersize = 10,
                 markerfacecolor = 'white',)
        scale = 0.5
        if (self.force_x != 0 or self.force_y != 0):
            force = math.hypot(self.force_x, self.force_y)
            force_x_norm = self.force_x / force
            force_y_norm = self.force_y / force
            plt.arrow(self.x,
                      self.y,
                      force_x_norm * scale,
                      force_y_norm * scale,
                      head_width = 0.15)
        
    def draw_offset(self, scale = 1):
        plt.plot(self.x + self.u * scale, 
                 self.y + self.v * scale, 
                 marker = 'o', 
                 markeredgecolor = 'k',
                 markersize = 10,
                 markerfacecolor = 'white',)

class Beam():
    """
    Балки
    """
    def __init__(self, point_0:Point, point_1:Point):
        self.points = [point_0, point_1]
        x_0 = self.points[0].x
        y_0 = self.points[0].y
        x_1 = self.points[1].x
        y_1 = self.points[1].y
        self.alpha = math.atan2(y_1 - y_0, x_1 - x_0)
    
    def draw(self):
        x_0 = self.points[0].x
        y_0 = self.points[0].y
        x_1 = self.points[1].x
        y_1 = self.points[1].y
        plt.plot([x_0, x_1], 
                 [y_0, y_1], 
                 linestyle = '-', 
                 
                 color = 'k', 
                 linewidth = 3)
        
    def draw_offset(self, scale = 1):
        x_0 = self.points[0].x + self.points[0].u * scale
        y_0 = self.points[0].y + self.points[0].v * scale
        x_1 = self.points[1].x + self.points[1].u * scale
        y_1 = self.points[1].y + self.points[1].v * scale
        plt.plot([x_0, x_1], 
                 [y_0, y_1], 
                 linestyle = '--', 
                 color = 'k', 
                 linewidth = 3)
        
    def length(self):
        x_0 = self.points[0].x
        y_0 = self.points[0].y
        x_1 = self.points[1].x
        y_1 = self.points[1].y
        return math.hypot(x_1 - x_0, y_1 - y_0)
        
    def define_stiffness_matrix(self, E, A):
        s = math.sin(self.alpha)
        c = math.cos(self.alpha)
        l = self.length()
        temp_sm = np.array(([[c * c, s * c, -c * c, -s * c],
                             [s * c, s * s, -s * c, -s * s],
                             [-c * c, -s * c, c * c, s * c],
                             [-s * c, - s * s, s * c, s * s]]))
        self.stiffness_matrix = temp_sm / l * E * A
        return
    
    def add_force(self, force_x, force_y, t):
        F = np.array([force_x, force_y])
        x_0 = self.points[0].x
        y_0 = self.points[0].y
        x_1 = self.points[1].x
        y_1 = self.points[1].y
        tau = np.array([x_1 - x_0, y_1 - y_0]) / math.hypot(x_1 - x_0, y_1 - y_0)
        F_tau = F * tau
        F_n = F - F_tau
        F_1tau = F_2tau = F_tau / 2
        F_1n = F_n * (1 - t)
        F_2n = F_n * t
        self.points[0].force_x = F_1n[0] + F_1tau[0]
        self.points[0].force_y = F_1n[1] + F_1tau[1]
        self.points[1].force_x = F_2n[0] + F_2tau[0]
        self.points[1].force_y = F_2n[1] + F_2tau[1]

class Bridge():
    """
    Мост
    """
    def __init__(self, list_of_points, list_of_beams):
        
        beams = []
        for beam in list_of_beams:
            point_0 = list_of_points[beam[0]]
            point_1 = list_of_points[beam[1]]
            beams.append(Beam(point_0, point_1))
        
        self.list_of_beams = list_of_beams
        self.points = list_of_points
        self.beams = beams
        
    def draw(self):
        for beam in self.beams:
            beam.draw()
        for point in self.points:
            point.draw()
            
    def draw_offset(self, scale = 1):
        for beam in self.beams:
            beam.draw_offset(scale)
        for point in self.points:
            point.draw_offset(scale)
    
    def make_glob_mk(self):
        E = 2.1e11
        A = 18 / 10000
        n = len(self.points)
        glob_mk = np.zeros((n*2, n*2))
        for i in range(len(self.beams)):
            self.beams[i].define_stiffness_matrix(E, A)
            el_mk = self.beams[i].stiffness_matrix
            m, k = self.list_of_beams[i]
            for u in range(2):
                for v in range(2):
                    glob_mk[m*2 + u, m*2 + v] += el_mk[u,v]
                    glob_mk[m*2 + u, k*2 + v] += el_mk[u,v+2]
                    glob_mk[k*2 + u, m*2 + v] += el_mk[u+2,v]
                    glob_mk[k*2 + u, k*2 + v] += el_mk[u+2,v+2]
                    
        for i, point in enumerate(self.points):
            if point.is_fixed:
                glob_mk[i*2, i*2] *= 1e12
                glob_mk[i*2 + 1, i*2 +1] *= 1e12
        self.glob_mk = glob_mk
        
    def make_force_vector(self):
        n = len(self.points)
        force_vector = np.zeros((n*2, 1))
        for i, point in enumerate(self.points):
            force_vector[i*2] = point.force_x
            force_vector[i*2 + 1] = point.force_y
        self.force_vector = force_vector
    
    def define_offset(self):
        delta = np.linalg.solve(self.glob_mk, self.force_vector)
        for i, point in enumerate(self.points):
            point.u = delta[i*2]
            point.v = delta[i*2 + 1]
            
    def clear_forces(self):
        for point in self.points:
            point.force_x = 0;
            point.force_y = 0;
        

            
    
def task_1():
    
    L = 12
    n = 7
    h = L/n*math.sqrt(3)/2
#    h = 4
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.set_ylim([-0.5, h + 0.5])
    points_lower_x = np.linspace(0, L, n+1)
    points_upper_x = np.linspace(0, L, n, endpoint = False) + L/n/2 
    list_of_points = [Point(points_lower_x[i], 0) for i in range(n+1)]
    list_of_points += [Point(points_upper_x[i], h) for i in range(n)]
#    for point in list_of_points:
#        point.draw()
        
    #Горизонтальные нижние
    l1 = [(i, i+1) for i in range(n)]
    #Горизонтальные верхние
    l2 = [(i+n+1, i+n+2) for i in range(n-1)]
    #диагонали /
    l3 = [(i, i+n+1) for i in range(n)]
    #диагонали \
    l4 = [(i+1, i+n+1) for i in range(n)]
    list_of_beams = l1 + l2 + l3 + l4
    
    bridge = Bridge(list_of_points, list_of_beams)
    
    bridge.draw()
    
def task_2():
    
    L = 12
    n = 10
    h = L/n*math.sqrt(3)/2
#    h = 4
    fig, ax = plt.subplots()
    ax.set_aspect(1)
    ax.set_ylim([-2, h + 0.5])
    points_lower_x = np.linspace(0, L, n+1)
    points_upper_x = np.linspace(0, L, n, endpoint = False) + L/n/2 
    list_of_points = [Point(points_lower_x[i], 0) for i in range(n+1)]
    list_of_points += [Point(points_upper_x[i], h) for i in range(n)]
    
    list_of_points[0].is_fixed = True
    list_of_points[n].is_fixed = True
    
    list_of_points[n // 2].force_y = -5000
    list_of_points[n // 2 +1].force_y = -5000
    list_of_points[n // 2 + 2].force_x = 0
    
#    for point in list_of_points:
#        point.draw()
        
    #Горизонтальные нижние
    l1 = [(i, i+1) for i in range(n)]
    #Горизонтальные верхние
    l2 = [(i+n+1, i+n+2) for i in range(n-1)]
    #диагонали /
    l3 = [(i, i+n+1) for i in range(n)]
    #диагонали \
    l4 = [(i+1, i+n+1) for i in range(n)]
    list_of_beams = l1 + l2 + l3 + l4
    
    list_of_points[n // 2].force_y = -5000
    list_of_points[n // 2 +1].force_y = -5000
    list_of_points[n // 2 + 2].force_x = 0
    
    bridge = Bridge(list_of_points, list_of_beams)
    
#    bridge.beams[3].add_force(10000, 10000, 0.5)
    bridge.draw()
    bridge.make_glob_mk()
    bridge.make_force_vector()
    bridge.define_offset()
    bridge.draw_offset(scale = 1000)
#    fig2, ax2 = plt.subplots()
#    plt.spy(bridge.glob_mk)
    
    print(np.shape(bridge.glob_mk))
#    print(np.linalg.matrix_rank(bridge.glob_mk))


def task_3():
    L = 12
    n = 10
    h = L/n*math.sqrt(3)/2

    points_lower_x = np.linspace(0, L, n+1)
    points_upper_x = np.linspace(0, L, n, endpoint = False) + L/n/2 
    list_of_points = [Point(points_lower_x[i], 0) for i in range(n+1)]
    list_of_points += [Point(points_upper_x[i], h) for i in range(n)]
    
    list_of_points[0].is_fixed = True
    list_of_points[n].is_fixed = True
        
    #Горизонтальные нижние
    l1 = [(i, i+1) for i in range(n)]
    #Горизонтальные верхние
    l2 = [(i+n+1, i+n+2) for i in range(n-1)]
    #диагонали /
    l3 = [(i, i+n+1) for i in range(n)]
    #диагонали \
    l4 = [(i+1, i+n+1) for i in range(n)]
    list_of_beams = l1 + l2 + l3 + l4
       
    bridge = Bridge(list_of_points, list_of_beams)
    
    L = 12
    n = 10
    h = L/n*math.sqrt(3)/2
#    h = 4
    
    points_lower_x = np.linspace(0, L, n+1)
    points_upper_x = np.linspace(0, L, n, endpoint = False) + L/n/2 
    list_of_points = [Point(points_lower_x[i], 0) for i in range(n+1)]
    list_of_points += [Point(points_upper_x[i], h) for i in range(n)]
    
    list_of_points[0].is_fixed = True
    list_of_points[n].is_fixed = True
    

        
    #Горизонтальные нижние
    l1 = [(i, i+1) for i in range(n)]
    #Горизонтальные верхние
    l2 = [(i+n+1, i+n+2) for i in range(n-1)]
    #диагонали /
    l3 = [(i, i+n+1) for i in range(n)]
    #диагонали \
    l4 = [(i+1, i+n+1) for i in range(n)]
    list_of_beams = l1 + l2 + l3 + l4

    bridge = Bridge(list_of_points, list_of_beams)
    
    a = L / n
    frames = 50
    mass = 500
    bridge.make_glob_mk()
    
    for i in range(frames):
        fig, ax = plt.subplots()
        ax.set_aspect(1)
        ax.set_ylim([-2, h + 0.5])
        
        bridge.clear_forces();
        b = L / frames * i
        number_of_beam = int(b // a)
        t = (b % a) / a
        bridge.beams[number_of_beam].add_force(0, - mass * 9.81, t)
        bridge.make_force_vector()
        bridge.define_offset()
        bridge.draw()
        bridge.draw_offset(scale = 1000)
        plt.savefig(fname = 'bridge' + str(i) + '.png')
        plt.close()

task_3()
    


        
        
