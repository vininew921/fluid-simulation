use rand::Rng;
use raylib::prelude::*;
use rayon::prelude::*;
use std::f32::consts::PI;

pub const GRAVITY: f32 = -180.0;
pub const COLLISION_DAMPING: f32 = 0.1;
pub const SMOOTHING_RADIUS: f32 = 10.0;
pub const MASS: f32 = 1.0;
pub const TARGET_DENSITY: f32 = 9.0;
pub const PRESSURE_MULTIPLIER: f32 = 90.0;
pub const PARTICLE_COUNT: usize = 600;
pub const PARTICLE_SIZE: f32 = 3.0;

pub struct ParticleSystem {
    pub positions: [Vector2; PARTICLE_COUNT],
    pub predicted_positions: [Vector2; PARTICLE_COUNT],
    pub velocities: [Vector2; PARTICLE_COUNT],
    pub densities: [f32; PARTICLE_COUNT],
    pub particle_properties: [f32; PARTICLE_COUNT],
    pub bounds: Vector2,
}

impl ParticleSystem {
    pub fn new(screen_width: i32, screen_height: i32, bounds: Vector2) -> Self {
        let mut positions = [Vector2::default(); PARTICLE_COUNT];
        let velocities = [Vector2::default(); PARTICLE_COUNT];
        let predicted_positions = [Vector2::default(); PARTICLE_COUNT];
        let densities = [0.0; PARTICLE_COUNT];
        let particle_properties = [0.0; PARTICLE_COUNT];

        for i in 0..PARTICLE_COUNT {
            positions[i] = Vector2 {
                x: raylib::get_random_value::<i32>(5, screen_width - 5) as f32,
                y: raylib::get_random_value::<i32>(5, screen_height - 5) as f32,
            };
        }

        ParticleSystem {
            positions,
            predicted_positions,
            velocities,
            densities,
            particle_properties,
            bounds,
        }
    }

    pub fn update(&mut self, delta_time: f32) {
        let mut predicted_positions = self.predicted_positions.clone();

        self.velocities
            .par_iter_mut()
            .zip(&mut predicted_positions)
            .enumerate()
            .for_each(|(i, (velocity, p_pos))| {
                *velocity += Vector2::new(0.0, -1.0) * GRAVITY * delta_time;
                *p_pos = self.positions[i] + *velocity * delta_time;
            });

        self.predicted_positions = predicted_positions.clone();

        for i in 0..PARTICLE_COUNT {
            self.densities[i] = self.calculate_density(self.predicted_positions[i]);
        }

        //Calculate and apply pressure forces
        for i in 0..PARTICLE_COUNT {
            let pressure_force = self.calculate_pressure_force(i);
            let pressure_acceleration = pressure_force / self.densities[i];
            self.velocities[i] += pressure_acceleration * delta_time;
        }

        //Update position and resolve collisions
        self.positions
            .par_iter_mut()
            .zip(&mut self.velocities)
            .for_each(|(position, velocity)| {
                *position += *velocity * delta_time;

                if position.x < 1.0 {
                    position.x = PARTICLE_SIZE + 0.2;
                    velocity.x *= -COLLISION_DAMPING;
                }
                if position.x >= self.bounds.x {
                    position.x = self.bounds.x - PARTICLE_SIZE;
                    velocity.x *= -COLLISION_DAMPING;
                }
                if position.y < 1.0 {
                    position.y = PARTICLE_SIZE + 0.2;
                    velocity.y *= -COLLISION_DAMPING;
                }
                if position.y >= self.bounds.y {
                    position.y = self.bounds.y - PARTICLE_SIZE;
                    velocity.y *= -COLLISION_DAMPING;
                }
            });
    }

    fn smoothing_kernel(&mut self, dst: f32, radius: f32) -> f32 {
        if dst >= radius {
            return 0.0;
        }

        let volume = (PI * f32::powf(radius, 4.0)) / 6.0;
        (radius - dst) * (radius - dst) / volume
    }

    fn smoothing_kernel_derivative(&mut self, dst: f32, radius: f32) -> f32 {
        if dst >= radius {
            return 0.0;
        }

        let scale = 12.0 / (PI * f32::powf(radius, 4.0));
        (dst - radius) * scale
    }

    fn calculate_density(&mut self, point: Vector2) -> f32 {
        let mut density = 0.0;

        for i in 0..PARTICLE_COUNT {
            let dst = (self.predicted_positions[i] - point).length();
            let influence = self.smoothing_kernel(dst, SMOOTHING_RADIUS);
            density += MASS * influence;
        }

        density
    }

    fn calculate_pressure_force(&mut self, particle_index: usize) -> Vector2 {
        let mut pressure_force = Vector2::zero();

        for other_particle_index in 0..PARTICLE_COUNT {
            if particle_index == other_particle_index {
                continue;
            }

            let offset = self.predicted_positions[other_particle_index]
                - self.predicted_positions[particle_index];

            let dst = offset.length();
            let dir = if dst == 0.0 {
                self.get_random_direction()
            } else {
                offset / dst
            };

            let slope = self.smoothing_kernel_derivative(dst, SMOOTHING_RADIUS);
            let density = self.densities[other_particle_index];
            let shared_pressure =
                self.calculate_shared_pressure(density, self.densities[particle_index]);
            pressure_force += -dir * shared_pressure * slope * MASS / density;
        }

        pressure_force
    }

    fn calculate_shared_pressure(&mut self, density_a: f32, density_b: f32) -> f32 {
        let pressure_a = self.convert_density_to_pressure(density_a);
        let pressure_b = self.convert_density_to_pressure(density_b);

        (pressure_a + pressure_b) / 2.0
    }

    fn convert_density_to_pressure(&mut self, density: f32) -> f32 {
        let density_error = density - TARGET_DENSITY;
        let pressure = PRESSURE_MULTIPLIER * density_error;
        pressure
    }

    fn get_random_direction(&mut self) -> Vector2 {
        let mut rng = rand::thread_rng();
        let dir = Vector2::new(rng.gen::<f32>() - 0.5, rng.gen::<f32>() - 0.5);
        dir.normalized()
    }
}
