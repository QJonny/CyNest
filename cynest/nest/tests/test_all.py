#! /usr/bin/env python
#
# test_all.py
#
# This file is part of NEST.
#
# Copyright (C) 2004 The NEST Initiative
#
# NEST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# NEST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with NEST.  If not, see <http://www.gnu.org/licenses/>.

import test_errors
import test_stack
import test_create
import test_status
import test_connectapi
import test_findconnections
import test_connectoptions
import test_events
import test_networks
import test_threads
import test_use
import test_dataconnect
import test_simulate

def run():
    test_errors.run()
    test_stack.run()
    test_create.run()                    
    test_status.run()
    test_connectapi.run()
    test_findconnections.run()    
    test_connectoptions.run()    
    test_events.run()
    test_networks.run()
    test_threads.run()  
    test_use.run()
    test_dataconnect.run()
    test_simulate.run()

