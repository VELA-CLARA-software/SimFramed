-- --------------------------------------------------------
-- Host:                         astecnas2
-- Server version:               5.5.62-MariaDB - Source distribution
-- Server OS:                    Linux
-- HeidiSQL Version:             10.2.0.5599
-- --------------------------------------------------------

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET NAMES utf8 */;
/*!50503 SET NAMES utf8mb4 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;


-- Dumping database structure for master_lattice
CREATE DATABASE IF NOT EXISTS `master_lattice` /*!40100 DEFAULT CHARACTER SET utf8 */;
USE `master_lattice`;

-- Dumping structure for table master_lattice.cavity
CREATE TABLE IF NOT EXISTS `cavity` (
  `name` varchar(128) NOT NULL,
  `cavity_type` text,
  `field_amplitude` double NOT NULL DEFAULT '0',
  `frequency` double NOT NULL DEFAULT '2998500000',
  `phase` double NOT NULL DEFAULT '0',
  `cell_length` double NOT NULL DEFAULT '0.033333',
  `field_definition` text,
  `field_definition_sdds` text,
  `field_definition_gdf` text,
  `lsc_enable` int(11) NOT NULL DEFAULT '1',
  `longitudinal_wakefield_sdds` text,
  `transverse_wakefield_sdds` text,
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.dipole
CREATE TABLE IF NOT EXISTS `dipole` (
  `name` varchar(128) NOT NULL,
  `angle` double NOT NULL DEFAULT '0',
  `length` double NOT NULL DEFAULT '0',
  `entrance_edge_angle` double NOT NULL DEFAULT '0',
  `exit_edge_angle` double NOT NULL DEFAULT '0',
  `edge_field_integral` double NOT NULL DEFAULT '0.4',
  `half_gap` double NOT NULL DEFAULT '0.02',
  `k1` double NOT NULL DEFAULT '0',
  `csr_bins` int(11) NOT NULL DEFAULT '20',
  `sr_enable` tinyint(4) NOT NULL DEFAULT '1',
  `isr_enable` double NOT NULL DEFAULT '1',
  `csr_output_filename` double DEFAULT NULL,
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 ROW_FORMAT=COMPACT;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.element
CREATE TABLE IF NOT EXISTS `element` (
  `name` varchar(128) NOT NULL DEFAULT '',
  `type` text NOT NULL,
  `parent` text NOT NULL,
  `length` float NOT NULL DEFAULT '0',
  `start_x` float NOT NULL DEFAULT '0',
  `start_y` float NOT NULL DEFAULT '0',
  `start_z` float NOT NULL DEFAULT '0',
  `end_x` float NOT NULL DEFAULT '0',
  `end_y` float NOT NULL DEFAULT '0',
  `end_z` float NOT NULL DEFAULT '0',
  `phi` float NOT NULL DEFAULT '0',
  `psi` float NOT NULL DEFAULT '0',
  `theta` float NOT NULL DEFAULT '0',
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.kicker
CREATE TABLE IF NOT EXISTS `kicker` (
  `name` varchar(128) NOT NULL,
  `horizontal_angle` double NOT NULL DEFAULT '0',
  `vertical_angle` double NOT NULL DEFAULT '0',
  `length` double NOT NULL DEFAULT '0',
  `sr_enable` int(11) NOT NULL DEFAULT '0',
  `isr_enable` int(11) NOT NULL DEFAULT '0',
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 ROW_FORMAT=COMPACT;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.laser
CREATE TABLE IF NOT EXISTS `laser` (
  `name` varchar(128) NOT NULL,
  `type` text,
  `value` double DEFAULT NULL,
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.longitudinal_wakefield
CREATE TABLE IF NOT EXISTS `longitudinal_wakefield` (
  `name` varchar(128) NOT NULL,
  `cell_length` double NOT NULL DEFAULT '0.033333',
  `field_definition` text,
  `scale_kick` double NOT NULL DEFAULT '1',
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 ROW_FORMAT=COMPACT;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.quadrupole
CREATE TABLE IF NOT EXISTS `quadrupole` (
  `name` varchar(128) NOT NULL,
  `k1l` double NOT NULL DEFAULT '0',
  `length` double NOT NULL DEFAULT '0',
  `sr_enable` int(11) NOT NULL DEFAULT '0',
  `isr_enable` int(11) NOT NULL DEFAULT '0',
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.screen
CREATE TABLE IF NOT EXISTS `screen` (
  `name` varchar(50) NOT NULL,
  `length` double NOT NULL DEFAULT '0',
  `output_filename` text,
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 ROW_FORMAT=COMPACT;

-- Data exporting was unselected.

-- Dumping structure for table master_lattice.solenoid
CREATE TABLE IF NOT EXISTS `solenoid` (
  `name` varchar(128) NOT NULL,
  `field_definition` text,
  `field_amplitude` double NOT NULL DEFAULT '0',
  `smooth` int(11) NOT NULL DEFAULT '0',
  UNIQUE KEY `name` (`name`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8 ROW_FORMAT=COMPACT;

-- Data exporting was unselected.

/*!40101 SET SQL_MODE=IFNULL(@OLD_SQL_MODE, '') */;
/*!40014 SET FOREIGN_KEY_CHECKS=IF(@OLD_FOREIGN_KEY_CHECKS IS NULL, 1, @OLD_FOREIGN_KEY_CHECKS) */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
