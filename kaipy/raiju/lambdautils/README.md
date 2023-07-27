
Set of scripts to generate lambda distributions for raiju
It seems like overkill (maybe it is), but the idea is to have a structured way of taking any previously generated raijuconfig.h5 file and know exactly how it was generated. That way, e.g. if we want to tweak the energy bounds or change the distribution between two values, we have a way of doing that while making sure everything else stays the same

The overall procedure for specifying a lambda distribution is:

1) Fill out as many specParams as you want (one for each species)
1a) First choose and fill out a DistType for each one, then create your specParams object
2) Give a list of your specParams to AlamData.py:AlamData. It will use the params to generate a list of Species objects, which will contain the fully realized lambda distributions for each alamParam
3) Call fileIO.SaveAlamConfig(...) to write the alam configuration to file