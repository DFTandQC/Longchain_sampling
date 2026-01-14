# Monomer Molecule Library

This folder stores monomer (monomer) molecule structure files used for cluster sampling.

## Current Files

- `opt-PT-B97-3c.xyz` - Platinum cluster monomer unit

## Adding New Molecules

To perform mixed cluster sampling containing other molecules, place the molecule's XYZ format structure file in this directory.

### Example: Adding H₂SO₄ Molecule

1. Optimize molecular structure using quantum chemistry software (e.g., ORCA, Gaussian, etc.)
2. Export to XYZ format: `h2so4.xyz`
3. Place in this folder
4. Reference in configuration file:
   ```json
   {
     "molecules": [
       {"name": "H2SO4", "file": "monomer/h2so4.xyz", "count": 1},
       {"name": "PT", "file": "monomer/opt-PT-B97-3c.xyz", "count": 2}
     ]
   }
   ```

## Important Notes

- XYZ file format: First line contains number of atoms, second line is comment/title, followed by atomic coordinates
- File paths can be absolute or relative to the project root directory
- **Molecule flexibility**: The system supports diverse molecule types:
  - Oxygen-containing molecules (acids, alcohols, ethers): O atoms define the head region
  - Nitrogen-containing molecules (amines, ammonia): N atoms define the head region
  - Sulfur-containing molecules (thiols, sulfides): S atoms define the head region
  - Hydrocarbons and other molecules: All atoms used for head region definition
- Examples of supported molecules:
  - ✅ H₂SO₄ (sulfuric acid) - Contains O atoms
  - ✅ NH₃ (ammonia) - Contains N atom
  - ✅ H₂O (water) - Contains O atom
  - ✅ Dimethyl ether - Contains O atom
  - ✅ Thiols - Contains S atom
- Verify that molecule files are not corrupted before use (test with: `head -3 filename.xyz`)
