Repository that chronicles the code involved at each stage of designing a plant to produce BioMethanol for my capstone project

The project is structured as a Julia package ecosystem to ensure separation between the physical calculations and plant data that must be hard coded.

In each sub-repo, there is a dedicated module for defining any constants that serves as the base of calculations. These usually are the lowest level dependecies for other more integrated modules. Therefore, any changes made in "hard coded" data will automatically update all calculations further downstream. 
