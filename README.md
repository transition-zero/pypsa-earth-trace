<!--
SPDX-FileCopyrightText:  PyPSA-Earth and PyPSA-Eur Authors

SPDX-License-Identifier: AGPL-3.0-or-later
-->

# PyPSA-Earth. A Flexible Python-based Open Optimisation Model to Study Energy System Futures around the World.

<p align="left">
by
<a href="https://pypsa-meets-earth.github.io">
    <img src="https://github.com/pypsa-meets-earth/pypsa-meets-earth.github.io/raw/main/assets/img/logo.png" width="150">
<a/>
</p>

## Development Status: **Stable and Active**

[![Status Linux](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-linux.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-linux.yaml)
[![Status Mac](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-mac.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-mac.yaml)
[![Status Windows](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-windows.yaml/badge.svg?branch=main&event=push)](https://github.com/pypsa-meets-earth/pypsa-earth/actions/workflows/ci-windows.yaml)
[![Documentation Status](https://readthedocs.org/projects/pypsa-earth/badge/?version=latest)](https://pypsa-earth.readthedocs.io/en/latest/?badge=latest)
![Size](https://img.shields.io/github/repo-size/pypsa-meets-earth/pypsa-earth)
[![License: AGPL v3](https://img.shields.io/badge/License-AGPLv3-blue.svg)](https://www.gnu.org/licenses/agpl-3.0)
[![REUSE status](https://api.reuse.software/badge/github.com/pypsa-meets-earth/pypsa-earth)](https://api.reuse.software/info/github.com/pypsa-meets-earth/pypsa-earth)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/pypsa-meets-earth/pypsa-earth/main.svg)](https://results.pre-commit.ci/latest/github/pypsa-meets-earth/pypsa-earth/main)
[![Discord](https://img.shields.io/discord/911692131440148490?logo=discord)](https://discord.gg/AnuJBk23FU)
[![Google Drive](https://img.shields.io/badge/Google%20Drive-4285F4?style=flat&logo=googledrive&logoColor=white)](https://drive.google.com/drive/folders/1U7fgktbxlaGzWxT2C0-Xv-_ffWCxAKZz)

**PyPSA-Earth is the first open-source global energy system model with data in high spatial and temporal resolution.** It enables large-scale collaboration by providing a tool that can model the world energy system or any subset of it. This work is derived from the European [PyPSA-Eur](https://pypsa-eur.readthedocs.io/en/latest/) model using new data and functions. It is suitable for operational as well as combined generation, storage and transmission expansion studies. The model provides two main features: (1) customizable data extraction and preparation scripts with global coverage and (2) a [PyPSA](https://pypsa.readthedocs.io/en/latest/) energy modelling framework integration. The data includes electricity demand, generation and medium to high-voltage networks from open sources, yet additional data can be further integrated. A broad range of clustering and grid meshing strategies help adapt the model to computational and practical needs.

The model is described in the Applied Energy article "PyPSA-Earth. A new global open energy system optimization model demonstrated in Africa", 2023, https://doi.org/10.1016/j.apenergy.2023.121096 [(BibTeX)](https://pypsa-earth.readthedocs.io/en/latest/talks_and_papers.html#publications). The [documentation](https://pypsa-earth.readthedocs.io/en/latest/index.html) provides additional information.

**PyPSA meets Earth is a free and open source software initiative aiming to develop a powerful energy system model for Earth.** We work on open data, open source modelling, open source solver support and open communities. Stay tuned and join our mission - We look for users, co-developers and leaders! Check out our [website for results and our projects](https://pypsa-meets-earth.github.io/projects.html). Happy coding!

<p align="center">
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/ddf041d1b98ca8f8c310f1c6393ec426ab5594cf.png" width=30%>
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/940b2673cfc31c4a6f01b7908f546d39d67df27e.png" width=23%>
  <img src="https://forum.openmod.org/uploads/db8804/original/1X/6af089c376b19b72ad148e4e4326c162b94db68f.png" width=35%>
</p>

<p align="center"><b> Figure:</b> Example power systems build with PyPSA-Earth. See images of ~193 more countries at <a>https://zenodo.org/records/10080766</a></p>

## Livetracker. Most popular global models:

<p align="center">
<a href="https://star-history.com/#pypsa-meets-earth/pypsa-earth&OSeMOSYS/osemosys_global&niclasmattsson/Supergrid&SGIModel/MUSE_OS&etsap-TIMES/TIMES_model&Date">
    <img src="https://api.star-history.com/svg?repos=pypsa-meets-earth/pypsa-earth,OSeMOSYS/osemosys_global,niclasmattsson/Supergrid,SGIModel/MUSE_OS,etsap-TIMES/TIMES_model&type=Date" width="60%">
<a/>

## Get involved

There are multiple ways to get involved and learn more about our work. That's how we organise ourselves:

- [**Discord NEW! (Open)**](https://discord.gg/AnuJBk23FU)
  - chat with the community, team up on features, exchange with developers, code in voice channels
  - registration and usage is for free
        <p align="left">
            <a href="https://discord.gg/AnuJBk23FU">
            <img src="https://discord.com/assets/cb48d2a8d4991281d7a6a95d2f58195e.svg" width="150">
            <a/>
        </p>
- **General initiative meeting (Open)**
  - every forth Thursday each month Thursday 16-17:00 (UK time)
    <a href="https://drive.google.com/file/d/1naH4WwW9drkOkOJ3PLO4fyWdkZQi5-_w/view?usp=share_link">
    `download .ics`
    </a>
  - join for project news and high-level code updates
  - meeting hosted on Discord
  - [open agenda](https://docs.google.com/document/d/1r6wm2RBe0DWFngmItpFfSFHA-CnUmVcVTkIKmthdW3g/edit?usp=sharing). See what we will discuss. Invited members have edit rights.
- **Buddy talk (Open)**
  - book a 30min meeting with Max to discuss anything you like
  - booking link: [calendly.com/pypsa-meets-earth](https://calendly.com/max-parzen/pypsa-meets-earth-exchange-30min)
- **Specific code meeting (Open)**
  - meeting hosted on Discord
  - join updates, demos, Q&A's, discussions and the coordination of each work package
    1. Demand creation and prediction meeting, on demand
    2. AI asset detection meeting, on demand
    3. Sector coupling meeting, every Thursday 09:00 (UK time), <a href="https://drive.google.com/file/d/1PDdmjsKhzyGRo0_YrP4wPQkn2XTNh6jA/view?usp=share_link" >`download .ics`</a>
    4. PyPSA-Earth meeting, every Thursday 16:00 (UK time), <a href="https://drive.google.com/file/d/1gaLmyV4qGPXsogkeRcAPWjC0ESebUxU-/view?usp=share_link" >`download .ics`</a>
- **Outreach meeting (Open)**
  - every second week, Tuesday 17:00 (UK time)
  - planning, discussing events, workshops, communication, community activities
- [**Google Drive**](https://drive.google.com/drive/folders/13Z8Y9zgsh5IZaDNkkRyo1wkoMgbdUxT5?usp=sharing)
  - access to minutes, presentations, lists, documents (access to minutes)

## Installation

1. Open your terminal at a location where you want to install pypsa-earth. Type the following in your terminal to download the package from GitHub:

   ```bash
      .../some/path/without/spaces % git clone https://github.com/pypsa-meets-earth/pypsa-earth.git
   ```
2. The python package requirements are curated in the `envs/environment.yaml` file.
   The environment can be installed using:

```bash
    .../pypsa-earth % conda env create -f envs/environment.yaml
```

   If the above takes longer than 30min, you might want to try mamba for faster installation:

```bash
    (base) conda install -c conda-forge mamba

    .../pypsa-earth % mamba env create -f envs/environment.yaml
```

3. For running the optimization one has to install the solver. We can recommend the open source HiGHs solver which installation manual is given [here](https://github.com/PyPSA/PyPSA/blob/633669d3f940ea256fb0a2313c7a499cbe0122a5/pypsa/linopt.py#L608-L632).
4. To use jupyter lab (new jupyter notebooks) **continue** with the [ipython kernel installation](http://echrislynch.com/2019/02/01/adding-an-environment-to-jupyter-notebooks/) and test if your jupyter lab works:

   ```bash
      .../pypsa-earth % ipython kernel install --user --name=pypsa-earth
      .../pypsa-earth % jupyter lab
   ```
5. Verify or install a java redistribution from the [official website](https://www.oracle.com/java/technologies/downloads/) or equivalent.
   To verify the successful installation the following code can be tested from bash:

   ```bash
      .../pypsa-earth % java -version
   ```

   The expected output should resemble the following:

   ```bash
      java version "1.8.0_341"
      Java(TM) SE Runtime Environment (build 1.8.0_341-b10)
      Java HotSpot(TM) 64-Bit Server VM (build 25.341-b10, mixed mode)
   ```

## Test run on tutorial

- In the folder open a terminal/command window to be located at this path `~/pypsa-earth/`
- Activate the environment `conda activate pypsa-earth`
- Rename config.tutorial.yaml to config.yaml. For instance in Linux:
  ```bash
  mv config.tutorial.yaml config.yaml
  ```
- Run a dryrun of the Snakemake workflow by typing simply in the terminal:
  ```bash
  snakemake -j 1 solve_all_networks -n
  ```

  Remove the -n to do a real run. Follow the tutorial of PyPSA-Eur 1 and 2 on [YouTube](https://www.youtube.com/watch?v=ty47YU1_eeQ) to continue with an analysis.

## Training

- We recently updated some [hackathon material](https://github.com/pypsa-meets-earth/documentation) for PyPSA-Earth. The hackathon contains jupyter notebooks with exercises. After going through the 1 day theoretical and practical material you should have a suitable coding setup and feel confident about contributing.
- The get a general feeling about the PyPSA functionality, we further recommend going through the [PyPSA](https://github.com/PyPSA/PyPSA/tree/master/examples) and [Atlite](https://github.com/PyPSA/atlite/tree/master/examples) examples.

## Questions and Issues

- We are happy to answer questions and help with issues **if they are public**. Through being public the wider community can benefit from the raised points. Some tips. **Bugs** and **feature requests** should be raised in the [**GitHub Issues**](https://github.com/pypsa-meets-earth/pypsa-earth/issues/new/choose). **General workflow** or **user questions** as well as discussion points should be posted at the [**GitHub Discussions**](https://github.com/pypsa-meets-earth/pypsa-earth/discussions/categories/q-a) tab. Happy coding.

## Documentation

The documentation is available here: [documentation](https://pypsa-earth.readthedocs.io/en/latest/index.html).

## Collaborators

<!-- https://github.com/marketplace/actions/contribute-list -->

<!-- readme: collaborators,contributors,restyled-commits/- -start -->
<table>
<tr>
    <td align="center">
        <a href="https://github.com/hazemful">
            <img src="https://avatars.githubusercontent.com/u/26235356?v=4" width="100;" alt="hazemful"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/fneum">
            <img src="https://avatars.githubusercontent.com/u/29101152?v=4" width="100;" alt="fneum"/>
            <br />
            <sub><b>Fabian Neumann</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/ekatef">
            <img src="https://avatars.githubusercontent.com/u/30229437?v=4" width="100;" alt="ekatef"/>
            <br />
            <sub><b>Ekaterina</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/euronion">
            <img src="https://avatars.githubusercontent.com/u/42553970?v=4" width="100;" alt="euronion"/>
            <br />
            <sub><b>Euronion</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Justus-coded">
            <img src="https://avatars.githubusercontent.com/u/44394641?v=4" width="100;" alt="Justus-coded"/>
            <br />
            <sub><b>Justus Ilemobayo</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/mnm-matin">
            <img src="https://avatars.githubusercontent.com/u/45293386?v=4" width="100;" alt="mnm-matin"/>
            <br />
            <sub><b>Mnm-matin</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/martacki">
            <img src="https://avatars.githubusercontent.com/u/53824825?v=4" width="100;" alt="martacki"/>
            <br />
            <sub><b>Martha Frysztacki</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/LukasFrankenQ">
            <img src="https://avatars.githubusercontent.com/u/55196140?v=4" width="100;" alt="LukasFrankenQ"/>
            <br />
            <sub><b>Lukas Franken</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/pz-max">
            <img src="https://avatars.githubusercontent.com/u/61968949?v=4" width="100;" alt="pz-max"/>
            <br />
            <sub><b>Max Parzen</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/davide-f">
            <img src="https://avatars.githubusercontent.com/u/67809479?v=4" width="100;" alt="davide-f"/>
            <br />
            <sub><b>Davide-f</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/koen-vg">
            <img src="https://avatars.githubusercontent.com/u/74298901?v=4" width="100;" alt="koen-vg"/>
            <br />
            <sub><b>Koen Van Greevenbroek</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/hazemakhalek">
            <img src="https://avatars.githubusercontent.com/u/87850910?v=4" width="100;" alt="hazemakhalek"/>
            <br />
            <sub><b>Hazem</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/energyLS">
            <img src="https://avatars.githubusercontent.com/u/89515385?v=4" width="100;" alt="energyLS"/>
            <br />
            <sub><b>EnergyLS</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/AnasAlgarei">
            <img src="https://avatars.githubusercontent.com/u/101210563?v=4" width="100;" alt="AnasAlgarei"/>
            <br />
            <sub><b>AnasAlgarei</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/yerbol-akhmetov">
            <img src="https://avatars.githubusercontent.com/u/113768325?v=4" width="100;" alt="yerbol-akhmetov"/>
            <br />
            <sub><b>Yerbol Akhmetov</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/DeniseGiub">
            <img src="https://avatars.githubusercontent.com/u/113139589?v=4" width="100;" alt="DeniseGiub"/>
            <br />
            <sub><b>DeniseGiub</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/GbotemiB">
            <img src="https://avatars.githubusercontent.com/u/48842684?v=4" width="100;" alt="GbotemiB"/>
            <br />
            <sub><b>Emmanuel Bolarinwa</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Tomkourou">
            <img src="https://avatars.githubusercontent.com/u/5240283?v=4" width="100;" alt="Tomkourou"/>
            <br />
            <sub><b>Thomas Kouroughli</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/GridGrapher">
            <img src="https://avatars.githubusercontent.com/u/127969728?v=4" width="100;" alt="GridGrapher"/>
            <br />
            <sub><b>GridGrapher</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Emre-Yorat89">
            <img src="https://avatars.githubusercontent.com/u/62134151?v=4" width="100;" alt="Emre-Yorat89"/>
            <br />
            <sub><b>Emre_Yorat</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/virio-andreyana">
            <img src="https://avatars.githubusercontent.com/u/114650479?v=4" width="100;" alt="virio-andreyana"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/giacfalk">
            <img src="https://avatars.githubusercontent.com/u/36954873?v=4" width="100;" alt="giacfalk"/>
            <br />
            <sub><b>Giacomo Falchetta</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Ekaterina-Vo">
            <img src="https://avatars.githubusercontent.com/u/99509555?v=4" width="100;" alt="Ekaterina-Vo"/>
            <br />
            <sub><b>Ekaterina-Vo</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/finozzifa">
            <img src="https://avatars.githubusercontent.com/u/167071962?v=4" width="100;" alt="finozzifa"/>
            <br />
            <sub><b>Finozzifa</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/cpschau">
            <img src="https://avatars.githubusercontent.com/u/124347782?v=4" width="100;" alt="cpschau"/>
            <br />
            <sub><b>Cschau</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Tooblippe">
            <img src="https://avatars.githubusercontent.com/u/805313?v=4" width="100;" alt="Tooblippe"/>
            <br />
            <sub><b>Tobias</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/doneachh">
            <img src="https://avatars.githubusercontent.com/u/132910766?v=4" width="100;" alt="doneachh"/>
            <br />
            <sub><b>Anton Achhammer</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/carlosfv92">
            <img src="https://avatars.githubusercontent.com/u/103258059?v=4" width="100;" alt="carlosfv92"/>
            <br />
            <sub><b>Carlos Fernandez</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/asolavi">
            <img src="https://avatars.githubusercontent.com/u/131155817?v=4" width="100;" alt="asolavi"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/stephenjlee">
            <img src="https://avatars.githubusercontent.com/u/11340470?v=4" width="100;" alt="stephenjlee"/>
            <br />
            <sub><b>Stephen J Lee</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/juli-a-ko">
            <img src="https://avatars.githubusercontent.com/u/126512394?v=4" width="100;" alt="juli-a-ko"/>
            <br />
            <sub><b>Juli-a-ko</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/squoilin">
            <img src="https://avatars.githubusercontent.com/u/4547840?v=4" width="100;" alt="squoilin"/>
            <br />
            <sub><b>Sylvain Quoilin</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/siddharth-krishna">
            <img src="https://avatars.githubusercontent.com/u/10712637?v=4" width="100;" alt="siddharth-krishna"/>
            <br />
            <sub><b>Siddharth Krishna</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/pitmonticone">
            <img src="https://avatars.githubusercontent.com/u/38562595?v=4" width="100;" alt="pitmonticone"/>
            <br />
            <sub><b>Pietro Monticone</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/Netotse">
            <img src="https://avatars.githubusercontent.com/u/89367243?v=4" width="100;" alt="Netotse"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/kma33">
            <img src="https://avatars.githubusercontent.com/u/25573938?v=4" width="100;" alt="kma33"/>
            <br />
            <sub><b>Katherine M. Antonio</b></sub>
        </a>
    </td></tr>
<tr>
    <td align="center">
        <a href="https://github.com/jessLryan">
            <img src="https://avatars.githubusercontent.com/u/122939887?v=4" width="100;" alt="jessLryan"/>
            <br />
            <sub><b>Jess</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/jarry7">
            <img src="https://avatars.githubusercontent.com/u/27745389?v=4" width="100;" alt="jarry7"/>
            <br />
            <sub><b>Jarrad Wright</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/HanaElattar">
            <img src="https://avatars.githubusercontent.com/u/87770004?v=4" width="100;" alt="HanaElattar"/>
            <br />
            <sub><b>HanaElattar</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/EmreYorat">
            <img src="https://avatars.githubusercontent.com/u/93644024?v=4" width="100;" alt="EmreYorat"/>
            <br />
            <sub><b>EmreYorat</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/AndreCNF">
            <img src="https://avatars.githubusercontent.com/u/19359510?v=4" width="100;" alt="AndreCNF"/>
            <br />
            <sub><b>André Cristóvão Neves Ferreira</b></sub>
        </a>
    </td>
    <td align="center">
        <a href="https://github.com/AlexanderMeisinger">
            <img src="https://avatars.githubusercontent.com/u/91368938?v=4" width="100;" alt="AlexanderMeisinger"/>
            <br />
            <sub><b>Null</b></sub>
        </a>
    </td></tr>
</table>
<!-- readme: collaborators,contributors,restyled-commits/- -end -->
