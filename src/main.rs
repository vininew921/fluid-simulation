pub mod particle_system;

use particle_system::ParticleSystem;
use raylib::prelude::*;

fn main() {
    let (screen_width, screen_height) = (400, 400);
    let bounds = Vector2::new(screen_width as f32, screen_height as f32);

    let mut particle_system = ParticleSystem::new(screen_width, screen_height, bounds);

    let (mut rl, thread) = raylib::init()
        .size(screen_width, screen_height)
        .title("Fluid simulation")
        .resizable()
        .build();

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::BLACK);
        d.draw_rectangle_lines(0, 0, bounds.x as i32, bounds.y as i32, Color::RED);

        particle_system.update(d.get_frame_time());

        //Draw particle
        for i in 0..particle_system::PARTICLE_COUNT {
            d.draw_circle_v(
                particle_system.positions[i],
                particle_system::PARTICLE_SIZE,
                Color::BLUE,
            );
        }
    }
}
