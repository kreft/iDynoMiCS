

A few general points about protocol files:
- There is NO automatic checking for consistency in the protocol file.
  Sometimes an inconsistency will cause an error in a simulation, and the
  code will attempt to tell you where in the protocol file the error
  occurred. However, the code WILL NOT catch errors in mass balances or
  units, so you must take care.

  Remember: "Garbage In, Garbage Out"

- Repeat: there is NO automatic checking for consistency - YOU must be careful!

- In general, for entering parameter values it is best to use the units given
  in the protocol file already, rather than change the described unit. The
  common units used are hour or day for time, g/L for concentrations, fg for
  mass, and micro-meter for spatial scale.

- Repeat: to be safe, CONVERT ALL PARAMETERS TO THE UNITS USED IN THIS EXAMPLE.

- Some keywords are reserved for special use by the software, such as:
  biomass, inert, capsule, pressure. 

- When you create a name for a variable, you are allowed to use upper and lower
  case and numbers, but no other symbols (and no spaces, either). Examples of
  valid names include: MyParameter, anotherParameter, parameter3.

- Any label that is prefixed with 'My' in the examples is free to be changed
  by you, the user, but take care not to use reserved keywords unless you
  are sure the usage is correct.


  Authors:
	Laurent Lardon
	Andreas Doetsch
	Jan-Ulrich Kreft
	Brian Merkey
	Cristian Picioreanu
	Barth F. Smets
	Joao Xavier
	Sonia Martins

  Website: http://www.idynomics.org
