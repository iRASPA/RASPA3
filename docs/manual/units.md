
# Units and conventions
\page units Units

## The standard units in `RASPA` from which all other units are derived are:

| Quantity      | Symbol   | Unit           | Value                                |
|---------------|----------|----------------|--------------------------------------|
| Length        | \f$l\f$      | Ångstrom        | \f$10^{-10}\f$ m                         |
| Temperature   | \f$T\f$      | Kelvin          | K                                    |
| Mass          | \f$m\f$      | Atomic mass     | \f$1.6605402 \times 10^{-27}\f$ kg       |
| Time          | \f$t\f$      | Pico seconds    | \f$10^{-12}\f$ s                         |
| Charge        | \f$q\f$      | Atomic charge   | \f$1.60217733 \times 10^{-19}\f$ C/particle |

## Some examples of derived units:

| Quantity             | Symbol   | Units                                                        | Conversion value                              |
|----------------------|----------|--------------------------------------------------------------|----------------------------------------------|
| Energy               | \f$U\f$      | \f$J = \text{mass} \times \text{length}^2 / \text{time}^2\f$     | \f$1.66054 \times 10^{-23}\f$ (=\f$10\f$ J/mol)      |
| Pressure             | \f$p\f$      | \f$\text{Pa} = \text{mass} / (\text{length} \times \text{time}^2)\f$ | \f$1.66054 \times 10^7\f$                        |
| Diffusion constant   | \f$D\f$      | \f$D = \text{length}^2 / \text{time}\f$                          | \f$1 \times 10^{-8}\f$                           |
| Force                | \f$f\f$      | \f$f = \text{length} / \text{time}^2\f$                          | \f$1.66054 \times 10^{-13}\f$                    |
| ...                  | ...      | ...                                                          | ...                                          |

A pressure input of 10 Pascal in the input file is converted to 'internal units' by dividing by \f$1.66054 \times 10^7\f$. In the output, any internal pressure is printed, multiplied by \f$1.66054 \times 10^7\f$. It is not necessary to convert units besides input and output, with a few exceptions.

One of the exceptions is the Coulombic conversion factor:
$$\frac{q_i q_j}{4\pi \epsilon_0}=\frac{\text{charge}^2}{4 \pi \times \text{electric constant} \times \text{length} \times \text{energy}} = 138935.4834964017$$
where the electric constant is \f$8.8541878176 \times 10^{-12}\f$ in units of \f$\text{C}^2/(\text{N} \cdot \text{m}^2)\f$. This factor is needed to convert the electrostatic energy to the internal units at every evaluation.

The Boltzmann constant \f$k_B\f$ is:
$$k_B = \text{Boltzmann constant}/\text{energy} = 0.8314464919$$
with the Boltzmann constant as \f$1.380650324 \times 10^{-23}\f$ in units of J/K, and \f$k_B = 0.8314464919\f$ in internal units.

-   Numbering is based on the `C`-convention, i.e. starting from zero.

-   Files in the current directory always have preference. Sometimes one would like to try various parameters for force field fitting, for example. In order to avoid making a lot of directories for each force field, it is more convenient to have the `force_field.def` file in the *current* directory.
