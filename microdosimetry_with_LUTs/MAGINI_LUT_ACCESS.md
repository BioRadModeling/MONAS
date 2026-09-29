Absolutely — here it is as **plain raw Markdown in one code block** so you can copy/paste it directly into `README.md`.

# Magini Lookup Tables

## 1. Request access

To request access to the Magini lookup tables, please send an email requesting access to the following individuals:

- **Giuseppe Schettino** — giuseppe.schettino@npl.co.uk
- **Francesco Romano** — francesco.romano@ct.infn.it
- **Dat Tran** — dtran2@lsu.edu
- **Marta Missiaggia** — mmissi2@lsu.edu

Please include the following information in your request:

- Your name
- Your institution or organization
- **Your GitHub username/account**
- A brief description of your intended use of the Magini LUTs

Your GitHub username is required because access to the Magini lookup tables is granted through a private GitHub repository.

Once your request has been reviewed and approved, the GitHub account provided in your request will be granted access to the private Magini LUT repository.

After access has been granted, continue with the installation instructions below to check out a local copy of the Magini lookup tables in the directory expected by MONAS.

## 2. Verify GitHub access

Before installing the lookup tables, make sure GitHub authentication is configured on your computer.

If you use SSH authentication, you can verify which GitHub account is being used with:

```bash
ssh -T git@github.com
```

A successful response will look similar to:

```text
Hi YOUR_GITHUB_USERNAME! You've successfully authenticated, but GitHub does not provide shell access.
```

Make sure that the GitHub username shown is the same account that was granted access to the private Magini LUT repository.

If the wrong GitHub account is displayed, update your SSH configuration or use the credentials associated with the GitHub account that was granted access.

## 3. Clone MONAS

If you have not already cloned MONAS:

```bash
git clone https://github.com/BioRadModeling/MONAS.git
cd MONAS
```

Check out the appropriate MONAS branch if necessary:

```bash
git switch feature_microdosimetry_from_LUTs
```

## 4. Install the Magini LUTs

The private Magini LUT repository is currently hosted at:

```text
github.com/dqtranPhysics/Magini-LUTs
```

The repository should be cloned directly into the location expected by MONAS.

From the root directory of your MONAS repository, run:

```bash
git clone git@github.com:dqtranPhysics/Magini-LUTs.git \
  microdosimetry_with_LUTs/lookup_tables/Magini
```

This command checks out the private Magini-LUTs repository directly as:

```text
microdosimetry_with_LUTs/lookup_tables/Magini/
```

The `.git/` directory inside `Magini/` belongs to the separate private Magini-LUTs repository.

The Magini LUT files are therefore not tracked as part of the public MONAS repository.

You can verify the installation with:

```bash
ls microdosimetry_with_LUTs/lookup_tables/Magini
```

You can also verify that the main lookup table exists with:

```bash
test -f microdosimetry_with_LUTs/lookup_tables/Magini/Magini.csv \
  && echo "Magini LUTs successfully installed"
```

If the installation was successful, you should see:

```text
Magini LUTs successfully installed
```

## 5. Access denied or repository not found

If you receive an error such as:

```text
ERROR: Repository not found.
fatal: Could not read from remote repository.
```

verify the following:

- Your request for Magini LUT access has been approved.
- The GitHub account you provided in your access request has been granted access to `dqtranPhysics/Magini-LUTs`.
- You are authenticated on your computer using that same GitHub account.
- The repository URL is entered correctly.

For SSH users, check which GitHub account is currently authenticated with:

```bash
ssh -T git@github.com
```

For example:

```text
Hi exampleUser! You've successfully authenticated, but GitHub does not provide shell access.
```

The username displayed by GitHub should match the GitHub account that was granted Magini LUT access.

If you have multiple GitHub accounts configured on the same computer, you may need to configure separate SSH keys or SSH host aliases so that Git uses the correct account.

## 6. Updating the Magini LUTs

The installed `Magini` directory is a separate Git repository.

If updated lookup tables are released, they can be retrieved without recloning MONAS.

From the MONAS root directory, run:

```bash
cd microdosimetry_with_LUTs/lookup_tables/Magini
git pull
```

This retrieves the latest authorized version of the Magini LUT repository.

You can then return to the MONAS root directory with:

```bash
cd ../../..
```

## 7. Removing and reinstalling the Magini LUTs

If you need to reinstall the Magini lookup tables, remove the local `Magini` directory:

```bash
rm -rf microdosimetry_with_LUTs/lookup_tables/Magini
```

Then clone the private repository again:

```bash
git clone git@github.com:dqtranPhysics/Magini-LUTs.git \
  microdosimetry_with_LUTs/lookup_tables/Magini
```

Be careful when using `rm -rf` and make sure you are removing only the intended `Magini` directory.

## Important

The Magini lookup tables are not part of the public MONAS repository and should not be committed to, redistributed through, or otherwise incorporated into the public MONAS source repository.

The public MONAS repository contains the software required to use the Magini method, while the lookup-table data themselves are distributed separately through the private:

```text
dqtranPhysics/Magini-LUTs
```

repository.

Access to the private Magini LUT repository is intended only for authorized users.

Please do not redistribute the lookup-table files or provide repository access to other users.

Individuals who wish to use the Magini LUTs should request access through the procedure described above.
