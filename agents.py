'''
Basic configuration of agents
'''

from enum import Enum
import uuid

class Status(Enum):
    '''
    Status of agents
    '''
    Susceptible = 's'
    Infected = 'i'
    Recovered = 'r'
    Death = 'd'

class Agent(object):
    '''
    agent class
    '''
    def __init__(self, x, y, velx, vely, size, status):
        self.id = int(uuid.uuid4())
        self.x = x
        self.y = y
        self.velx = velx
        self.vely = vely
        self.size = size
        self.status = status
        self.infected_time = 0

    