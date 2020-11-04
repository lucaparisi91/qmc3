{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 57,
   "metadata": {},
   "outputs": [],
   "source": [
    "import pandas as pd\n",
    "import json\n",
    "import os\n",
    "import numpy as np"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 58,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "(50, 3, 11)"
      ]
     },
     "execution_count": 58,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "dirname=\"testConfig\"\n",
    "with open(os.path.join(dirname,\"description.dat\")) as f:\n",
    "    settings=json.load(f)\n",
    "\n",
    "particleFile=os.path.join(dirname,\"particles.dat\")\n",
    "data=np.fromfile(particleFile,dtype=np.float64)\n",
    "data=data.reshape(settings[\"timeSlices\"],settings[\"dimensions\"],settings[\"nChains\"])\n",
    "data.shape"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 59,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "{'dimensions': 3,\n",
       " 'groups': [{'iEnd': 10, 'iEndExtended': 10, 'iStart': 0}],\n",
       " 'nChains': 11,\n",
       " 'timeSlices': 50,\n",
       " 'worms': [{'iChainHead': 10, 'iChainTail': 5, 'iHead': 10, 'iTail': 10}]}"
      ]
     },
     "execution_count": 59,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "settings"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 52,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "(1650,)"
      ]
     },
     "execution_count": 52,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "data.shape"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 48,
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1650.0"
      ]
     },
     "execution_count": 48,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "13200/8"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.8.5"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 4
}
